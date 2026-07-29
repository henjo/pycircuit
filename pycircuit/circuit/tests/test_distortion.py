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


## ============ higher-order plan, stage A: the scalar series ============

"""Stage A of ``doc/distortion_higher_order_plan.md``.

Builds the perturbation terms genuinely order by order via the Faà di Bruno
composition formula -- the thing an earlier attempt substituted Picard
iteration for. Scalar and circuit-free by design, so the combinatorics can be
gated on their own before harmonics are involved.

The gate, declared in the plan before this code: for ``Y x + b x**2 = u``
every term must be a *single monomial* with Catalan coefficients
1, -1, 2, -5, 14. Anything that expands to a sum means the composition
bookkeeping is wrong.
"""

from pycircuit.circuit.distortion import (_compositions,
                                          composition_coefficient,
                                          perturbation_terms)


def _quadratic_series(order):
    """Terms of Y x + b x**2 = u, symbolically."""
    Y, b, u = sympy.symbols('Y b u')
    G = 1 / Y
    x0 = G * u
    derivatives = {1: 2 * b * x0, 2: 2 * b, 3: 0, 4: 0, 5: 0, 6: 0}
    return perturbation_terms(G, b * x0 ** 2, derivatives, x0, order)


def test_compositions_enumerates_ordered_tuples():
    """Ordered, not unordered -- collapsing them would drop multiplicities."""
    assert sorted(_compositions(3, 1)) == [(3,)]
    assert sorted(_compositions(3, 2)) == [(1, 2), (2, 1)]
    assert sorted(_compositions(3, 3)) == [(1, 1, 1)]
    assert sorted(_compositions(4, 2)) == [(1, 3), (2, 2), (3, 1)]
    for k in range(1, 7):
        total = sum(len(list(_compositions(k, m))) for m in range(1, k + 1))
        assert total == 2 ** (k - 1), (
            'compositions of %d should number 2**(k-1)' % k)


def test_composition_formula_matches_a_direct_series_expansion():
    """Check the formula against sympy expanding f(x0+y) itself.

    Deliberately with a *concrete* f: a first version of this check used
    symbolic Derivative placeholders, whose substitution silently failed and
    reported the formula wrong from k=3. The formula was right; the check was
    not.
    """
    eps, x0, b, c = sympy.symbols('eps x0 b c')
    xs = sympy.symbols('x1:6')
    f = lambda v: b * v ** 2 + c * v ** 3
    derivatives = {m: sympy.diff(f(x0), x0, m) for m in range(1, 6)}
    terms = [x0] + list(xs)

    y = sum(eps ** n * xs[n - 1] for n in range(1, 6))
    direct = sympy.expand(sympy.series(f(x0 + y), eps, 0, 6).removeO())

    for k in range(1, 6):
        formula = composition_coefficient(k, terms, derivatives)
        assert sympy.simplify(sympy.expand(direct.coeff(eps, k) - formula)) == 0, (
            'composition formula disagrees at order %d' % k)


def test_scalar_terms_are_single_monomials_with_catalan_coefficients():
    """The stage A gate, exactly as declared in the plan."""
    Y, b, u = sympy.symbols('Y b u')
    terms = _quadratic_series(5)
    expected = [u / Y,
                -b * u ** 2 / Y ** 3,
                2 * b ** 2 * u ** 3 / Y ** 5,
                -5 * b ** 3 * u ** 4 / Y ** 7,
                14 * b ** 4 * u ** 5 / Y ** 9,
                -42 * b ** 5 * u ** 6 / Y ** 11]

    for order, (got, want) in enumerate(zip(terms, expected)):
        assert len(sympy.Add.make_args(sympy.expand(got))) == 1, (
            'x(%d) is not a single monomial: %s' % (order, sympy.expand(got)))
        assert sympy.simplify(got - want) == 0, (
            'x(%d): got %s, expected %s' % (order, sympy.simplify(got), want))


def test_partial_sums_converge_to_the_exact_root():
    """Numerically, the truncations must approach the true solution.

    ``Y x + b x**2 = u`` has a closed-form root, so this needs no simulation
    and no reference implementation.
    """
    Y, b, u = sympy.symbols('Y b u')
    values = {Y: 2.0, b: 0.7, u: 1.0}
    terms = [float(t.subs(values)) for t in _quadratic_series(4)]
    exact = (-2.0 + np.sqrt(4.0 + 4 * 0.7 * 1.0)) / (2 * 0.7)

    errors = [abs(sum(terms[:k + 1]) - exact) for k in range(5)]
    assert all(errors[i] > errors[i + 1] for i in range(4)), (
        'partial sums should improve monotonically, got %s'
        % ['%.3e' % e for e in errors])


def test_the_series_is_not_the_picard_iteration():
    """Guard the distinction this whole plan exists because of.

    An earlier attempt reported Picard results as perturbation results. These
    must not silently become the same object again.
    """
    Y, b, u = sympy.symbols('Y b u')
    values = {Y: 2.0, b: 0.7, u: 1.0}
    terms = [float(t.subs(values)) for t in _quadratic_series(3)]

    G, B, U = 1 / 2.0, 0.7, 1.0
    picard = G * U
    for _ in range(2):
        picard = G * (U - B * picard ** 2)

    assert abs(picard - sum(terms[:3])) > 1e-6, (
        'the second Picard iterate and the second-order truncation must '
        'differ -- they are different constructions')


## ====== higher-order plan, stage B: harmonics + graded truncation ======

"""Stage B of ``doc/distortion_higher_order_plan.md``.

Adds the harmonic index and grading by power of the drive, then truncates at
``U**3``. The gate: the result must reproduce ``harmonic_response`` exactly,
because the published second-order form *is* the consistent ``U**3``
truncation. If it does not, the grading is wrong and no higher-order number
built on it can be trusted.

It does, on harmonics 2 and 3, to floating point. The fundamental differs by
design -- ``harmonic_response`` drops its ``U**3`` correction at assembly (see
"The fundamental is not corrected" in the theory page) while the graded form
carries it -- and that difference is verified to scale as ``U**2`` relative to
the fundamental, which is exactly the signature of an ``O(U**3)`` term.
"""

from pycircuit.circuit.distortion import GradedSpectrum, graded_response

_GR_Y0, _GR_C0, _GR_B, _GR_C = 4.4e-3, 1e-9, 0.35, 4.7
_GR_W0 = 2 * np.pi * 1e4


def _graded_setup(drive):
    Y = lambda s: _GR_Y0 + s * _GR_C0
    response = lambda s: 1.0 / Y(s)

    def apply_G(harmonic, rhs):
        return [np.asarray(rhs, dtype=complex)[0]
                / Y(1j * harmonic.frequency((_GR_W0,)))]

    published = harmonic_response(
        apply_G, [drive],
        CubicNonlinearity(np.array([[_GR_B]]), np.array([[_GR_C]])),
        tones=(_GR_W0,), toolkit=numeric)
    graded = graded_response(response, drive, [_GR_B, _GR_C], (_GR_W0,),
                             max_power=3)
    return published, graded


def test_graded_spectrum_products_match_the_analytic_coefficients():
    """The convolution must give textbook harmonic coefficients.

    For x = Re[X e^{jt}]: x**2 has second-harmonic coefficient X**2/2 and DC
    |X|**2/2; x**3 has third-harmonic X**3/4 and fundamental 3|X|**2 X/4.
    """
    x = GradedSpectrum.from_phasor(1, 1, 2 + 0j)
    square, cube = x * x, x * x * x

    assert abs(square.phasor(2) - 2.0) < 1e-15      # X**2/2
    assert abs(square.phasor(0) - 2.0) < 1e-15      # |X|**2/2
    assert abs(cube.phasor(3) - 2.0) < 1e-15        # X**3/4
    assert abs(cube.phasor(1) - 6.0) < 1e-15        # 3|X|**2 X/4


def test_grading_tracks_the_power_of_the_drive():
    """Harmonic m from a single tone must first appear at exactly U**m."""
    x = GradedSpectrum.from_phasor(1, 1, 1 + 0j)
    assert x.powers_present(1) == [1]
    assert (x * x).powers_present(2) == [2]
    assert (x * x * x).powers_present(3) == [3]


@pytest.mark.parametrize('drive', [1e-5, 1e-6, 1e-7])
def test_graded_truncation_reproduces_the_published_form(drive):
    """The stage B gate, on the harmonics both constructions report."""
    published, graded = _graded_setup(drive)

    for m in (2, 3):
        want = published.amplitude((m,), 0)
        got = graded.phasor(m)
        assert abs(got - want) / abs(want) < 1e-12, (
            'harmonic %d at drive %g: graded %r vs published %r'
            % (m, drive, got, want))


def test_the_graded_form_carries_the_fundamental_correction():
    """And the published one does not -- verified as an O(U**3) difference.

    Not a disagreement: ``harmonic_response`` drops this term deliberately.
    Checking that the difference scales as U**2 relative to the fundamental
    confirms it is the O(U**3) correction and not an error somewhere.
    """
    ratios = []
    for drive in (1e-5, 1e-6):
        _, graded = _graded_setup(drive)
        linear = graded.phasor(1, max_power=1)
        ratios.append(abs(graded.phasor(1) - linear) / abs(linear))

    assert ratios[0] > 10 * ratios[1], (
        'the correction should fall as U**2, i.e. ~100x for a decade of '
        'drive; got %.3e -> %.3e' % tuple(ratios))


def test_the_test_point_is_inside_the_validity_range():
    """Guard the guard: these parameters must be weakly nonlinear.

    A first attempt at this comparison used a drive giving |X2/X1| = 0.90 --
    essentially at the divergence threshold -- where the 'correction' to the
    fundamental came out five times the fundamental. The numbers were
    arithmetically right and the operating point meaningless.
    """
    published, _ = _graded_setup(1e-5)
    assert published.HD2(0) < 0.15, (
        'test drive is not weakly nonlinear: HD2 = %.3f' % published.HD2(0))


## ====== higher-order plan, stage C: does raising the order help? ======

"""Stages C and D of ``doc/distortion_higher_order_plan.md``.

**It helps, substantially.** Against the transient cross-check on the biased
diode, raising the truncation from ``U**3`` to ``U**5`` improves accuracy at
every drive, and by two orders of magnitude in the middle of the valid range.

This is the opposite of what an earlier experiment concluded. That experiment
raised *Picard iterations* rather than the perturbation order and reported
non-monotonic behaviour; the genuine truncation is monotone. The tests below
pin the corrected result so the earlier one cannot quietly return.
"""

import math as _math


def _diode_taylor(highest_power):
    """Taylor coefficients of the biased diode, quadratic term upward."""
    VT = numeric.kboltzmann * 300 / numeric.qelectron
    v0 = brentq_operating_point()
    IS_eff = 1e-13 * np.exp(v0 / VT)
    return [IS_eff / (_math.factorial(n) * VT ** n)
            for n in range(2, highest_power + 1)]


def brentq_operating_point():
    from scipy.optimize import brentq
    VT = numeric.kboltzmann * 300 / numeric.qelectron
    return brentq(lambda v: v / 1e3 + 1e-13 * np.expm1(v / VT) - 1e-3, 0, 1)


def _diode_graded_hd(drive, max_power):
    VT = numeric.kboltzmann * 300 / numeric.qelectron
    v0 = brentq_operating_point()
    alpha = (1e-13 / VT) * np.exp(v0 / VT)
    w0 = 2 * np.pi * 1e4
    response = lambda s: 1.0 / (1 / 1e3 + s * 1e-9 + alpha)

    g = graded_response(response, drive, _diode_taylor(max_power), (w0,),
                        max_power=max_power)
    fundamental = abs(g.phasor(1))
    return abs(g.phasor(2)) / fundamental, abs(g.phasor(3)) / fundamental


@pytest.mark.parametrize('drive', [1e-4, 2e-4])
def test_higher_truncation_improves_accuracy(drive):
    """The stage C gate: U**5 must beat U**3 against the transient.

    Uses the measurement helper from test_distortion_vs_transient, so the
    reference is a real simulation and not another perturbation calculation.
    """
    from pycircuit.circuit.tests.test_distortion_vs_transient import _measure

    hd2_meas, hd3_meas, _, _ = _measure(drive)
    err = {}
    for power in (3, 5):
        hd2, hd3 = _diode_graded_hd(drive, power)
        err[power] = (abs(hd2_meas / hd2 - 1), abs(hd3_meas / hd3 - 1))

    assert err[5][0] < err[3][0], (
        'HD2 at drive %g: U**5 error %.4f should beat U**3 error %.4f'
        % (drive, err[5][0], err[3][0]))
    assert err[5][1] < err[3][1], (
        'HD3 at drive %g: U**5 error %.4f should beat U**3 error %.4f'
        % (drive, err[5][1], err[3][1]))


def test_the_improvement_is_large_enough_to_be_worth_it():
    """Stage D's gate: better than 2x inside the validity range.

    Declared in the plan before the measurement, so that "it helps a little"
    would have been recorded as a negative result rather than argued up.
    """
    from pycircuit.circuit.tests.test_distortion_vs_transient import _measure

    drive = 1e-4                      # perturbation ratio ~0.05, comfortably valid
    hd2_meas, hd3_meas, _, _ = _measure(drive)

    hd2_3, hd3_3 = _diode_graded_hd(drive, 3)
    hd2_5, hd3_5 = _diode_graded_hd(drive, 5)

    gain_hd2 = abs(hd2_meas / hd2_3 - 1) / max(abs(hd2_meas / hd2_5 - 1), 1e-12)
    assert gain_hd2 > 2.0, (
        'stage D wants better than 2x; HD2 improved only %.2fx' % gain_hd2)


def test_symbolic_growth_stays_polynomial_with_order():
    """The other half of stage C: the method must remain symbolically usable.

    Terms in the third harmonic go 2 -> 7 -> 16 for truncations U**3, U**5,
    U**7 -- polynomial in the order. An iterative scheme doubles instead (see
    the limits page), which is what would make it unusable symbolically.
    """
    Y0, b, c, d, e, A = sympy.symbols('Y0 b c d e A')
    w0 = sympy.Symbol('w0', positive=True)
    response = lambda s: 1 / Y0

    counts = []
    for power, coeffs in ((3, [b, c]), (5, [b, c, d, e])):
        g = graded_response(response, A, coeffs, (w0,), max_power=power)
        counts.append(len(sympy.Add.make_args(sympy.expand(g.phasor(3)))))

    assert counts[0] < counts[1] < 4 * counts[0], (
        'third-harmonic term count should grow polynomially, got %s' % counts)


# ---------------------------------------------------------------------------
# Multi-node graded response.  Gates declared in doc/distortion_mimo_plan.md.
# ---------------------------------------------------------------------------

import math
from scipy.optimize import brentq
from pycircuit.circuit.distortion import GradedVector, graded_response_mimo

## The 2013 paper's worked example 2: Tow-Thomas gm-C biquad, eq. (46).
## Centre 10.6931 MHz, Q = 20, unity gain at centre.
_BQ_G1 = 31.26e-6
## g_2 is NEGATIVE in the paper (p. 492: "g_2 = -31.26 uA/V"), verified at
## 600 dpi.  An earlier version of this file had it positive, which made the
## pencil come out with right-half-plane poles and produced a "defect in the
## published circuit" that did not exist.  The sign is load-bearing: it is
## the damping term of a Q = 20 resonator.
_BQ_G2 = -31.26e-6
_BQ_G3, _BQ_G4 = 625.2e-6, -625.2e-6
_BQ_C1 = _BQ_C2 = 9.3054e-12
_BQ_ALPHA = -0.0535                    # every cubic tied to its linear
_BQ_G1C, _BQ_G2C, _BQ_G3C, _BQ_G4C = (
    _BQ_ALPHA * g for g in (_BQ_G1, _BQ_G2, _BQ_G3, _BQ_G4))


def _bq_q(w):
    return _BQ_G3*_BQ_G4 + _BQ_C1*_BQ_C2*w**2 + _BQ_G2*_BQ_C2*1j*w


def _bq_solve(s, rhs):
    M = np.array([[-_BQ_G2 + s*_BQ_C1, -_BQ_G4],
                  [-_BQ_G3, s*_BQ_C2]], dtype=complex)
    return np.linalg.solve(M, np.asarray(rhs, dtype=complex))


def _cube(sp, maxp):
    return ((sp * sp).truncated(maxp) * sp).truncated(maxp)


def _bq_run(Xin, w, maxp, with_g4c=True, with_g1c=True):
    src = GradedSpectrum.from_phasor(1, 1, _BQ_G1*Xin)
    if with_g1c:
        src = src + _cube(GradedSpectrum.from_phasor(1, 1, Xin),
                          maxp).scaled(_BQ_G1C)

    def f(x):
        f1 = _cube(x[0], maxp).scaled(-_BQ_G2C)
        if with_g4c:
            f1 = f1 + _cube(x[1], maxp).scaled(-_BQ_G4C)
        return GradedVector([f1, _cube(x[0], maxp).scaled(-_BQ_G3C)])

    return graded_response_mimo(_bq_solve,
                                GradedVector([src, GradedSpectrum()]),
                                f, (w,), max_power=maxp)


def _bq_eq48(w, Xin, full):
    """The paper's published third harmonic, eq. (48).

    ``full=False`` keeps only the second fraction -- the ``g_3c``/``g_2c``
    part, which is what a hand elimination of node 2 can reach.  The first
    fraction carries ``g_4c`` and ``g_1c`` and is unreachable scalar-wise.
    """
    second = -Xin**3*_BQ_C2*1j*w*(
        _BQ_G1**3*_BQ_C2**2*w**2*(_BQ_G3C*_BQ_G4 + 3*_BQ_G2C*_BQ_C2*1j*w)
        ) / (4*_bq_q(3*w)*_bq_q(w)**3)
    if not full:
        return second
    first = 3*(_BQ_G4C*_BQ_G1**3*_BQ_G3**3 - _BQ_G1C*_bq_q(w)**3)
    return Xin**3*_BQ_C2*1j*w*first/(4*_bq_q(3*w)*_bq_q(w)**3) + second


def _scalar_poly(spec, coeffs, maxp):
    f = GradedSpectrum()
    power = spec
    for c in coeffs:
        power = (power * spec).truncated(maxp)
        if c != 0:
            f = f + power.scaled(c)
    return f


def test_mimo_contains_scalar_exactly():
    """Stage A gate: a 1x1 system must reproduce ``graded_response``.

    Not an approximation -- the multi-node path has to *contain* the scalar
    one, or no result computed through it can be trusted.  Checked across
    drives, truncation orders and harmonics rather than at one point, because
    the grading is what is under test and it is indexed by both.
    """
    IS_D = 1e-13
    VT = numeric.kboltzmann*300/numeric.qelectron
    Ib, Rl, Cl, f0 = 1e-3, 1e3, 1e-9, 1e4
    v0 = brentq(lambda v: v/Rl + IS_D*np.expm1(v/VT) - Ib, 0, 1)
    ISe = IS_D*np.exp(v0/VT)
    w0 = 2*np.pi*f0
    resp = lambda s: 1.0/(1/Rl + s*Cl + ISe/VT)
    taylor = [ISe/(math.factorial(n)*VT**n) for n in range(2, 14)]

    for A in (1e-4, 4e-4, 6e-4):
        for power in (3, 5, 7, 9, 11):
            coeffs = taylor[:power-1]
            want = graded_response(resp, A, coeffs, (w0,), max_power=power)
            got = graded_response_mimo(
                lambda s, rhs: [resp(s)*rhs[0]],
                GradedVector([GradedSpectrum.from_phasor(1, 1, A)]),
                lambda x: GradedVector([_scalar_poly(x[0], coeffs, power)]),
                (w0,), max_power=power)
            for m in (1, 2, 3, 5):
                a, b = want.phasor(m), got[0].phasor(m)
                assert abs(a - b) <= 1e-13*abs(a) + 1e-300, (
                    'MIMO != scalar at A=%g U^%d harmonic %d' % (A, power, m))


def test_mimo_carries_independent_problems_side_by_side():
    """Stage A gate, second half: a block-diagonal system must not mix nodes.

    A 1x1 check cannot catch a driver that accidentally sums right-hand sides
    across nodes before solving, because there is only one.
    """
    VT = numeric.kboltzmann*300/numeric.qelectron
    w0 = 2*np.pi*1e4
    r1 = lambda s: 1.0/(1e-3 + s*1e-9 + 0.016)
    r2 = lambda s: 1.0/(5e-4 + s*1e-9 + 0.016)
    c1, c2 = [0.3/VT**2, 0.2/VT**3], [0.15/VT**2, 0.1/VT**3]

    want1 = graded_response(r1, 3e-4, c1, (w0,), max_power=5)
    want2 = graded_response(r2, 5e-4, c2, (w0,), max_power=5)
    got = graded_response_mimo(
        lambda s, rhs: [r1(s)*rhs[0], r2(s)*rhs[1]],
        GradedVector([GradedSpectrum.from_phasor(1, 1, 3e-4),
                      GradedSpectrum.from_phasor(1, 1, 5e-4)]),
        lambda x: GradedVector([_scalar_poly(x[0], c1, 5),
                                _scalar_poly(x[1], c2, 5)]),
        (w0,), max_power=5)

    for i, want in enumerate((want1, want2)):
        for m in (1, 2, 3):
            a, b = want.phasor(m), got[i].phasor(m)
            assert abs(a - b) <= 1e-13*abs(a) + 1e-300


@pytest.mark.parametrize('freq', [1e5, 1e6, 5e6, 1.06931e7, 2e7])
def test_biquad_reachable_part_matches_published_eq48(freq):
    """Stage B gate: the matrix path must reproduce the scalar-reachable part.

    Isolates the matrix plumbing from the new nonlinearities: with ``g_4c``
    and ``g_1c`` switched off, the answer is one a hand elimination of node 2
    already reproduces, so any discrepancy is in the plumbing.
    """
    w = 2*np.pi*freq
    got = _bq_run(1e-3, w, 3, with_g4c=False, with_g1c=False)[0].phasor(3)
    want = _bq_eq48(w, 1e-3, full=False)
    assert abs(got - want) <= 1e-12*abs(want)


@pytest.mark.parametrize('freq', [1e5, 1e6, 5e6, 1.06931e7, 2e7])
def test_biquad_full_matches_published_eq48(freq):
    """Stage C gate: the complete published closed form, all four cubics.

    ``g_4c`` acts on node 2 and ``g_1c`` on the input, and neither is
    reachable by eliminating node 2 by hand: ``x2 = L(s) x1`` with ``L``
    frequency-dependent, so ``(L x1)**3 != L**3 x1**3`` once ``x1`` carries
    more than one harmonic.  This is new capability rather than a regression,
    and the gate is a formula the authors printed.
    """
    w = 2*np.pi*freq
    got = _bq_run(1e-3, w, 3)[0].phasor(3)
    want = _bq_eq48(w, 1e-3, full=True)
    assert abs(got - want) <= 1e-12*abs(want)


def test_biquad_omitted_nonlinearities_are_not_negligible():
    """Guards the *point* of the multi-node work, not just its arithmetic.

    If the terms only reachable through the matrix path were small, stage C
    would be a curiosity.  At 100 kHz they dominate by four orders of
    magnitude, so a scalar reduction is not a mild approximation there -- it
    misses essentially the whole answer.  A refactor that quietly dropped
    them would still pass a loose comparison; this fails.
    """
    w = 2*np.pi*1e5
    partial = abs(_bq_run(1e-3, w, 3, with_g4c=False,
                          with_g1c=False)[0].phasor(3))
    full = abs(_bq_run(1e-3, w, 3)[0].phasor(3))
    assert full > 1e3*partial, (
        'expected the node-2 and input cubics to dominate at low frequency, '
        'got full=%g partial=%g' % (full, partial))


## Stage D: the 2013 paper's worked example 1, 3-stage RNMC amplifier.
## Structurally different from the biquad: three nodes, *quadratic* as well
## as cubic terms, and nonlinearities that act on a node other than the one
## they are injected into.
_AMP = dict(gm1=245e-6, g01=1/98e3, gm2=200e-6, g02=1/107e3,
            gm2q=4e-3, gm2c=60e-3, gm3=1e-3, g03=1/23e3,
            gm3q=2e-3, gm3c=3e-3, CL=10e-12, C1=2e-12, C2=1e-12)


def _amp_p(w):
    return _AMP['g03'] + 1j*w*_AMP['CL']


def _amp_q(w):
    a = _AMP
    return (a['g01']*a['g02']*a['g03']
            + 1j*w*(a['CL']*a['g01']*a['g02']
                    + a['C1']*(a['g02']*a['g03'] + a['gm2']*a['gm3'])
                    + a['C2']*a['g03']*(a['g02'] + a['gm2']))
            - w**2*((a['C1']+a['C2'])*a['CL']*a['g02']
                    + a['C2']*a['CL']*a['gm2']))


def _amp_r(w):
    a = _AMP
    return a['g01']*a['g02'] + 1j*w*(a['C1']*a['g02']
                                     + a['C2']*(a['g02'] + a['gm2']))


def _amp_solve(s, rhs):
    """Open-loop form of published eq. (45).

    The paper prints the *feedback* configuration; the open-loop matrix is
    its limit ``k_in -> 1``, ``k_3 -> 0``, ``g_03e -> g_03``.  That
    reconstruction is not assumed -- it is checked against published eq. (39)
    by :func:`test_amplifier_matrix_reproduces_published_linear_solution`.
    """
    a = _AMP
    M = np.array([[a['g01'] + (a['C1']+a['C2'])*s, -a['C2']*s, -a['C1']*s],
                  [a['gm2'], a['g02'], 0.0],
                  [0.0, a['gm3'], -a['g03'] - a['CL']*s]], dtype=complex)
    return np.linalg.solve(M, np.asarray(rhs, dtype=complex))


def _amp_run(Xin, w, maxp):
    ## Published eqs. (41)/(42) are explicitly restricted to the output
    ## transconductor's nonlinearities, so gm2q/gm2c and g03q/g03c stay off.
    def f(x):
        sq = (x[1] * x[1]).truncated(maxp)
        cu = (sq * x[1]).truncated(maxp)
        return GradedVector([GradedSpectrum(), GradedSpectrum(),
                             sq.scaled(_AMP['gm3q']) + cu.scaled(_AMP['gm3c'])])

    src = GradedVector([GradedSpectrum.from_phasor(1, 1, _AMP['gm1']*Xin),
                        GradedSpectrum(), GradedSpectrum()])
    return graded_response_mimo(_amp_solve, src, f, (w,), max_power=maxp)


def _amp_linear_node3(w, Xin):
    """Published eq. (39), third node -- the linearised fundamental."""
    return -_AMP['gm1']*_AMP['gm2']*_AMP['gm3']*Xin/_amp_q(w)


def _amp_published(w, Xin):
    a = _AMP
    hd2 = abs(a['gm3q']*a['gm1']*a['gm2']*Xin*_amp_p(w)**2*_amp_r(2*w)
              / (2*a['gm3']*_amp_q(w)*_amp_q(2*w)))
    hd3 = (abs(a['gm1']**2*a['gm2']**2*Xin**2*_amp_p(w)**3*_amp_r(3*w)
               / (a['gm3']*_amp_q(w)**2*_amp_q(3*w)))
           * abs(a['gm3c']/4
                 - 1j*w*a['C1']*a['gm3q']**2*a['gm2']/_amp_q(2*w)))
    return hd2, hd3


@pytest.mark.parametrize('freq', [1e2, 1e4, 1e6, 1e7])
def test_amplifier_matrix_reproduces_published_linear_solution(freq):
    """The open-loop matrix is reconstructed, so it gets its own gate.

    The paper prints eq. (45) for the feedback configuration only; the
    open-loop matrix used here is its limit.  Published eq. (39) gives the
    linearised solution at all three nodes independently, so it checks the
    reconstruction without reference to any nonlinearity.
    """
    w = 2*np.pi*freq
    Xin = 1e-4
    got = _amp_solve(1j*w, [_AMP['gm1']*Xin, 0, 0])
    want = [_AMP['gm1']*_AMP['g02']*Xin*_amp_p(w)/_amp_q(w),
            -_AMP['gm1']*_AMP['gm2']*Xin*_amp_p(w)/_amp_q(w),
            _amp_linear_node3(w, Xin)]
    for node, (a, b) in enumerate(zip(got, want)):
        assert abs(a - b) <= 1e-12*abs(b), 'node %d' % (node + 1)


@pytest.mark.parametrize('freq', [1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
def test_amplifier_hd2_hd3_match_published_eq41_eq42(freq):
    """Stage D gate: a second circuit, and the first with quadratic terms.

    HD3 is the sign-sensitive one: it sums a direct cubic contribution
    against a second-order contribution of the *quadratic* coefficient, so a
    wrong sign convention fails it while HD2, a magnitude ratio, would still
    pass.  That is why both are checked and not just the larger.

    Ratios are taken against the **linearised** fundamental, which is what
    the paper's ``X_{1,3}`` denotes.  The graded fundamental additionally
    carries its own ``U**3`` correction; see the companion test.
    """
    w = 2*np.pi*freq
    Xin = 1e-4
    g = _amp_run(Xin, w, 3)
    lin = _amp_linear_node3(w, Xin)
    want2, want3 = _amp_published(w, Xin)
    assert abs(abs(g[2].phasor(2)/lin) - want2) <= 1e-11*want2
    assert abs(abs(g[2].phasor(3)/lin) - want3) <= 1e-11*want3


def test_amplifier_fundamental_correction_scales_as_drive_squared():
    """Pins *why* the published ratio uses the linearised fundamental.

    Dividing by the graded fundamental instead disagrees with eq. (41) by
    5.7e-3 at 100 Hz, which looks like an error and is not: the graded form
    carries the fundamental's ``U**3`` correction and the published one drops
    it.  The signature of that -- rather than of a mistake -- is that the gap
    falls by 100x per decade of drive.  Asserting the *scaling* rather than
    the value is what distinguishes the two explanations.
    """
    w = 2*np.pi*1e2
    gaps = []
    for Xin in (1e-4, 1e-5, 1e-6):
        g = _amp_run(Xin, w, 3)
        hd2 = abs(g[2].phasor(2)/g[2].phasor(1))
        want = _amp_published(w, Xin)[0]
        gaps.append(abs(hd2 - want)/want)

    for coarse, fine in zip(gaps, gaps[1:]):
        assert 90.0 < coarse/fine < 110.0, (
            'gap should scale as U**2, got ratios %s' % gaps)


# ---------------------------------------------------------------------------
# Stage E: higher order on a multi-node circuit.
# ---------------------------------------------------------------------------

def _bq_run_with(solve, Xin, w, maxp):
    src = (GradedSpectrum.from_phasor(1, 1, _BQ_G1*Xin)
           + _cube(GradedSpectrum.from_phasor(1, 1, Xin), maxp).scaled(_BQ_G1C))

    def f(x):
        return GradedVector([_cube(x[0], maxp).scaled(-_BQ_G2C)
                             + _cube(x[1], maxp).scaled(-_BQ_G4C),
                             _cube(x[0], maxp).scaled(-_BQ_G3C)])

    return graded_response_mimo(solve, GradedVector([src, GradedSpectrum()]),
                                f, (w,), max_power=maxp)


def _bq_centre():
    return np.sqrt(-_BQ_G3*_BQ_G4/(_BQ_C1*_BQ_C2))


def test_biquad_pencil_as_published_is_stable():
    """The published circuit is stable, and this pins the sign that makes it so.

    **This test replaces one that asserted the opposite.**  An earlier version
    of this file took ``g_2 = +31.26 uA/V``; the paper gives ``-31.26 uA/V``
    (p. 492, verified at 600 dpi).  With the wrong sign the pencil has
    right-half-plane poles, and that artifact was written up as a defect in
    the published circuit.  It was a defect in the transcription.

    Retained as a test because the sign is easy to lose again -- it is the
    damping term of a ``Q = 20`` resonator, and every magnitude comparison in
    this file passes either way, since those compare our arithmetic against
    the paper's formula evaluated with *the same* constants.  Nothing except
    a stability check or a time-domain integration notices.
    """
    assert _BQ_G2 < 0, 'g_2 is negative in the source'
    roots = np.roots([_BQ_C1*_BQ_C2, -_BQ_G2*_BQ_C2, -_BQ_G3*_BQ_G4])
    assert all(r.real < 0 for r in roots), roots

    ## Centre frequency and Q are unaffected by the sign, which is part of
    ## why it survived: the circuit looks entirely reasonable either way.
    w0 = _bq_centre()
    assert abs(w0/(2*np.pi)/1e6 - 10.6931) < 1e-3
    assert abs(w0*_BQ_C1/abs(_BQ_G2) - 20.0) < 1e-2


def test_magnitude_gates_cannot_see_the_damping_sign():
    """Why the wrong sign survived four gates at 1e-14.

    Every published comparison here evaluates the paper's formula with the
    same constants our code uses, so a wrong constant cancels out and the
    gate passes regardless.  That is worth a test of its own, because it
    bounds what those gates establish: they check the *machinery* against the
    *algebra*, and say nothing about whether the constants are the paper's.
    """
    w = _bq_centre()
    ours = abs(_bq_run(1e-3, w, 3)[0].phasor(3))
    theirs = abs(_bq_eq48(w, 1e-3, full=True))
    assert abs(ours - theirs) <= 1e-12*theirs

    ## Flip the damping sign in *both* and the agreement is untouched.
    global _BQ_G2
    saved = _BQ_G2
    try:
        _BQ_G2 = -saved
        flipped_ours = abs(_bq_run(1e-3, w, 3)[0].phasor(3))
        flipped_theirs = abs(_bq_eq48(w, 1e-3, full=True))
        assert abs(flipped_ours - flipped_theirs) <= 1e-12*flipped_theirs
    finally:
        _BQ_G2 = saved

    ## ...even though the circuit is a different one.
    assert abs(flipped_ours - ours) > 1e-6*ours


def test_multinode_symbolic_growth_stays_polynomial():
    """Stage E gate: the claim that this method suits a symbolic tool has to
    survive going to several nodes.

    Both sides use an ``s``-dependent response and a pure cubic, so the only
    difference is one node against two coupled ones.  A first attempt gave
    the scalar case ``response = 1/Y`` -- no ``s`` at all, so it could never
    accumulate denominators -- and measured raw ``count_ops``; that pair of
    mistakes made multi-node look like it blew up by 16x per two orders.  It
    does not.
    """
    a, c1, c2, b, d, k3, A, w0 = sympy.symbols('a c1 c2 b d k3 A w0')

    def solve(s, rhs):
        det = (a + s*c1)*(s*c2) - b*d
        return [(s*c2*rhs[0] + b*rhs[1])/det,
                (d*rhs[0] + (a + s*c1)*rhs[1])/det]

    def size(e):
        num, _ = sympy.fraction(sympy.together(e))
        return len(sympy.Add.make_args(sympy.expand(num)))

    ratios = []
    for power in (3, 5, 7):
        coeffs = [0, k3] + [0]*(power - 2)
        scalar = graded_response(lambda s: 1/(a + s*c1), A, coeffs,
                                 (w0,), max_power=power)
        mimo = graded_response_mimo(
            solve,
            GradedVector([GradedSpectrum.from_phasor(1, 1, A),
                          GradedSpectrum()]),
            lambda x: GradedVector([GradedSpectrum(),
                                    _cube(x[0], power).scaled(k3)]),
            (w0,), max_power=power)

        ## Keys are the structural measure: exactly one set per node.
        assert (sum(len(c.terms) for c in mimo.components)
                == 2*len(scalar.terms))
        ratios.append(size(mimo[0].phasor(3)) / size(scalar.phasor(3)))

    assert max(ratios) < 5.0, (
        'multi-node should cost a small constant factor, got %s' % ratios)


@pytest.mark.parametrize('Xin', [0.05, 0.10])
def test_biquad_converges_monotonically_below_the_bound(Xin):
    """Stage E gate: error must fall with every added order, below the bound.

    Reference is a direct integration of the biquad's own ODEs -- genuinely
    independent of the perturbation machinery, sharing only the circuit.  It
    is insensitive to cycles (200/400/800) and tolerance (1e-11..1e-13) to
    nine significant figures, so it is not the limiting factor here.

    At ``Xin = 0.3`` this ordering breaks (the sequence goes 1.0e-1, 1.3e-1,
    8.2e-3, 1.2e-2, 1.6e-3).  That is the convergence bound, not a defect,
    and it is why the drives here are below it.
    """
    from scipy.integrate import solve_ivp

    w = _bq_centre()

    def rhs(t, y):
        x1, x2 = y
        u = Xin*np.cos(w*t)
        return [(_BQ_G2*x1 + _BQ_G4*x2 + _BQ_G1*u + _BQ_G1C*u**3
                 + _BQ_G2C*x1**3 + _BQ_G4C*x2**3)/_BQ_C1,
                (_BQ_G3*x1 + _BQ_G3C*x1**3)/_BQ_C2]

    period = 2*np.pi/w
    tend = 200*period
    keep = 32
    ts = np.linspace(tend - keep*period, tend, keep*256, endpoint=False)
    sol = solve_ivp(rhs, (0, tend), [0.0, 0.0], t_eval=ts, method='DOP853',
                    rtol=1e-11, atol=1e-18, max_step=period/40)
    spectrum = np.fft.rfft(sol.y[0])/len(sol.y[0])
    want = 2*abs(spectrum[3*keep])

    errors = [abs(abs(_bq_run_with(_bq_solve, Xin, w, n)[0].phasor(3))
                  - want)/want
              for n in (3, 5, 7, 9, 11)]

    for coarse, fine in zip(errors, errors[1:]):
        assert fine < coarse, 'not monotone: %s' % errors
    assert errors[-1] < 1e-6, 'expected deep convergence, got %s' % errors


# ---------------------------------------------------------------------------
# Two tones on the multi-node path (plan section 8.1).
# ---------------------------------------------------------------------------

def _amp_run_two_tone(Xin, w1, w2, maxp=3):
    src = (GradedSpectrum.from_phasor((1, 0), 1, _AMP['gm1']*Xin)
           + GradedSpectrum.from_phasor((0, 1), 1, _AMP['gm1']*Xin))

    def f(x):
        sq = (x[1] * x[1]).truncated(maxp)
        cu = (sq * x[1]).truncated(maxp)
        return GradedVector([GradedSpectrum(), GradedSpectrum(),
                             sq.scaled(_AMP['gm3q']) + cu.scaled(_AMP['gm3c'])])

    return graded_response_mimo(
        _amp_solve,
        GradedVector([src, GradedSpectrum(), GradedSpectrum()]),
        f, (w1, w2), max_power=maxp)


def test_harmonic_index_accepts_ints_and_tuples_alike():
    """The tuple refactor must not change the single-tone API.

    Indices are stored as tuples internally so that two-tone convolution is
    componentwise addition, but ``phasor(3)`` has to keep meaning the third
    harmonic of one tone.  A bare int is normalised on the way in.
    """
    s = GradedSpectrum.from_phasor(2, 1, 4.0)
    assert s.phasor(2) == s.phasor((2,))
    assert s.powers_present(2) == s.powers_present((2,)) == [1]
    assert list(s.terms) == [((2,), 1), ((-2,), 1)]


def test_two_tone_index_arithmetic_is_componentwise():
    """``(2,-1)`` must arise as a plain sum of one-sided components.

    This is what storing harmonics two-sided buys: no case analysis over
    which combinations of sum and difference frequencies produce a given
    product.  A concatenating ``+`` (the default for tuples) would silently
    produce nonsense indices instead of failing.
    """
    a = GradedSpectrum.from_phasor((1, 0), 1, 2.0)
    b = GradedSpectrum.from_phasor((0, 1), 1, 2.0)
    cube = ((a + b) * (a + b) * (a + b))
    assert ((2, -1), 3) in cube.terms
    assert ((3, 0), 3) in cube.terms
    for index, _ in cube.terms:
        assert len(index) == 2, 'index %r has been concatenated' % (index,)


@pytest.mark.parametrize('f1', [1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
def test_amplifier_im3_matches_published_eq43(f1):
    """Plan section 8.1 gate: two tones, three nodes, against eq. (43).

    The tone ratio ``f2 = 0.9 f1`` and drive are the paper's own (its
    Fig. 4b).  As with eqs. (41)/(42), the formula is restricted to the
    output transconductor's nonlinearities and the ratio is taken against
    the *linearised* fundamental.

    The tildes in the printed formula (``p~(jw2)``, ``q~(jw2)``) are complex
    conjugates, because the ``(2,-1)`` product takes ``w2`` with a negative
    sign.  **They make no difference to this comparison** -- measured, after
    an earlier version of this docstring asserted they did: the factor is a
    pure product and ratio inside ``|.|``, and conjugating any factor leaves
    a modulus unchanged.  They are kept because the formula is transcribed
    as printed, not because the test can see them.

    That is the same reason the right-half-plane pencil of eq. (46) is
    invisible in this paper's results, and it is worth noticing twice: these
    closed forms are moduli of products, so they are blind to an entire class
    of sign and phase error.  Agreement with them is strong evidence about
    magnitudes and no evidence at all about phase.

    Six frequencies rather than one because they span six decades of IM3
    magnitude and the circuit's poles lie inside that range.
    """
    a = _AMP
    Xin = 1e-4
    w1 = 2*np.pi*f1
    w2 = 0.9*w1

    got = (abs(_amp_run_two_tone(Xin, w1, w2)[2].phasor((2, -1)))
           / abs(_amp_linear_node3(w1, Xin)))

    first = abs(a['gm1']**2*a['gm2']**2*Xin**2*_amp_p(w1)**2
                * np.conj(_amp_p(w2))*_amp_r(2*w1 - w2)
                / (a['gm3']*_amp_q(w1)*np.conj(_amp_q(w2))*_amp_q(2*w1 - w2)))
    second = abs(3*a['gm3c']/4
                 - a['C1']*a['gm3q']**2*a['gm2']
                 * (1j*w1/_amp_q(2*w1) + 1j*(w1 - w2)/_amp_q(w1 - w2)))
    want = first*second

    assert abs(got - want) <= 1e-11*want


# ---------------------------------------------------------------------------
# Phase.  Every published gate in this module compares moduli; this does not.
# ---------------------------------------------------------------------------

def test_biquad_phase_and_waveform_match_the_ode_integration():
    """Validates what the published closed forms structurally cannot.

    Every gate against Buonomo & Lo Schiavo compares ``|.|`` of a product and
    ratio of transfer functions.  Conjugating any factor leaves such a
    modulus unchanged, so those comparisons are blind to sign and phase --
    demonstrated twice in this file: the ``p~``/``q~`` conjugates in eq. (43)
    make no difference to the result, and the eq. (46) pencil has
    right-half-plane poles that eqs. (47)/(48) do not notice.

    A time-domain integration has no such blindness, because a waveform is
    not a modulus.  Three checks here, in increasing strength:

    1. per-harmonic phase against the integration,
    2. mutations that preserve ``|X|`` exactly, to show the phase criterion
       can actually fail,
    3. the reconstructed waveform against the integrated one, which tests
       every harmonic's magnitude *and* phase at once.

    The FFT window starts at an exact multiple of the drive period, so its
    phase reference is ``cos(w0 t)`` -- the same one the phasors use.  Get
    that wrong and every phase is offset by a constant, which is why check 2
    matters: it distinguishes "phases agree" from "phases agree up to an
    offset nobody checked".
    """
    from scipy.integrate import solve_ivp

    Xin = 0.10
    w = _bq_centre()

    def rhs(t, y):
        x1, x2 = y
        u = Xin*np.cos(w*t)
        return [(_BQ_G2*x1 + _BQ_G4*x2 + _BQ_G1*u + _BQ_G1C*u**3
                 + _BQ_G2C*x1**3 + _BQ_G4C*x2**3)/_BQ_C1,
                (_BQ_G3*x1 + _BQ_G3C*x1**3)/_BQ_C2]

    period = 2*np.pi/w
    tend = 200*period
    keep = 32
    ts = np.linspace(tend - keep*period, tend, keep*256, endpoint=False)
    sol = solve_ivp(rhs, (0, tend), [0.0, 0.0], t_eval=ts, method='DOP853',
                    rtol=1e-12, atol=1e-18, max_step=period/60)
    spectrum = np.fft.rfft(sol.y[0])/len(sol.y[0])
    want = {m: 2*spectrum[m*keep] for m in (1, 3, 5)}

    got = _bq_run_with(_bq_solve, Xin, w, 11)

    ## 1 -- phase, harmonic by harmonic.
    for m in (1, 3, 5):
        offset = np.angle(got[0].phasor(m)/want[m], deg=True)
        assert abs(offset) < 0.01, 'harmonic %d off by %.4f deg' % (m, offset)

    ## 2 -- the criterion must be able to fail.  Each of these leaves the
    ## magnitude bit-identical, so a magnitude-only gate accepts all three.
    reference = got[0].phasor(3)
    for mutate, expected in ((lambda z: -z, 180.0),
                             (lambda z: 1j*z, 90.0)):
        mutated = mutate(reference)
        assert abs(mutated) == abs(reference), 'mutation changed |X|'
        offset = abs(np.angle(mutated/want[3], deg=True))
        assert abs(offset - expected) < 0.1, (
            'expected %g deg, got %g' % (expected, offset))

    ## 3 -- the whole waveform, which folds in every harmonic at once.
    rebuilt = np.full_like(ts, float(np.real(got[0].phasor(0))))
    for m in range(1, 12):
        phasor = got[0].phasor(m)
        if phasor != 0:
            rebuilt = rebuilt + np.real(phasor*np.exp(1j*m*w*ts))
    peak = np.max(np.abs(sol.y[0]))
    assert np.max(np.abs(rebuilt - sol.y[0]))/peak < 1e-8
