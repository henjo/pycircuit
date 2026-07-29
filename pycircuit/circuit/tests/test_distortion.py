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
from pycircuit.circuit.distortion import (Harmonic, harmonic_response,
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
    sol = harmonic_response(_single_node_apply_G(Y), [A], bm, cm,
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

    sol = harmonic_response(apply_G, U, b, c, tones=(w,), toolkit=numeric)
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
