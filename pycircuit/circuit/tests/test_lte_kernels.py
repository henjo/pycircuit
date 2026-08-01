"""Tests for the shared integrator algebra (pycircuit.circuit._lte_kernels).

STAGE 9(a).  These formulas were transcribed between `integrator.py` and
`jaxtransient.py` rather than shared, and every defect stage 9 found was a repair
that existed on one copy and not the other -- the Gear-2 error constant (4i, found
again by gate 9-1(a)), the opening-step ramp (stage 3, found again by 9(g)), and
the breakpoint scan (9(d), wrong in two different ways at once).

The point of the module is that ONE definition is reached by both backends, so
these tests check that property rather than re-deriving the arithmetic: the
independent derivations live in `test_transient_repairs.py`, which computes the
exact truncation error from an analytic charge and must NOT be pointed at the
kernel or it would be testing the code against itself.
"""
import numpy as np
import pytest

from pycircuit.circuit import _lte_kernels as K


def test_kernel_module_imports_no_backend():
    """The kernels must be plain arithmetic, so they trace under jax.jit.

    A `numpy` call in here would work on the CPU path and silently break tracing
    on the JAX one -- the failure would appear as a JAX error in a file that never
    mentions JAX.  Asserted structurally rather than by trying to trace, so it
    fails for the right reason on a machine with no JAX installed.
    """
    import inspect
    src = inspect.getsource(K)
    for banned in ('import numpy', 'import jax', 'from numpy', 'from jax'):
        assert banned not in src, \
            '_lte_kernels must not import a backend, found %r' % banned


def test_bdf2_alphas_annihilate_a_quadratic():
    """The defining property, not the transcribed expression.

    The VSS Gear-2 coefficients differentiate the quadratic through the three
    points exactly, so applied to a quadratic charge they must return its exact
    derivative.  Checked on a NON-uniform grid, since fixed-step SPICE2 formulas
    pass a uniform-grid check and are wrong the moment the step changes -- which
    is the whole reason these are computed rather than tabulated.
    """
    h1, h2 = 1e-5, 3.7e-6          # deliberately unequal
    q = lambda t: 3.0 + 2.0 * t + 5.0 * t ** 2
    dq = lambda t: 2.0 + 10.0 * t
    t_n = 0.371e-6
    got = K.bdf2_derivative(q(t_n), q(t_n - h1), q(t_n - h1 - h2), h1, h2)
    assert got == pytest.approx(dq(t_n), rel=1e-9)


def test_companions_match_their_textbook_forms():
    """Each companion current, against its definition written out separately."""
    h, h2 = 1e-5, 7e-6
    q_n, q_1, q_2 = 1.7, 1.3, 0.8
    iq_1, C = 0.5, 2e-9

    iq, geq = K.euler_companion(q_n, C, q_1, h)
    assert iq == pytest.approx((q_n - q_1) / h, rel=1e-15)
    assert geq == pytest.approx(C / h, rel=1e-15)

    iq, geq = K.trapezoidal_companion(q_n, C, q_1, iq_1, h)
    assert iq == pytest.approx(2.0 * (q_n - q_1) / h - iq_1, rel=1e-15)
    assert geq == pytest.approx(2.0 * C / h, rel=1e-15)

    a0, a1, a2 = K.bdf2_alphas(h, h2)
    iq, geq = K.bdf2_companion(q_n, C, q_1, q_2, h, h2)
    assert iq == pytest.approx(a0 * q_n + a1 * q_1 + a2 * q_2, rel=1e-15)
    assert geq == pytest.approx(a0 * C, rel=1e-15)


def test_gear2_error_constant_is_one_third_not_one_quarter():
    """The 3/4 optimism, pinned in the one place the constant now lives.

    YWR's Table I GEAR2 residual gives (1/4) h^2 q''' against a true (1/3).  The
    CPU carried that until stage 4i and the JAX backend until gate 9-1(a) -- the
    same defect found three times in two transcriptions, which is why the constant
    is derived once here.
    """
    h, q3 = 1e-6, 1.0e6
    assert K.gear2_lte(h, h, q3 / 6.0) == pytest.approx(-(1.0 / 3.0) * h ** 2 * q3,
                                                        rel=1e-12)
    assert K.trapezoidal_lte(h, q3 / 6.0) == pytest.approx(-(1.0 / 6.0) * h ** 2 * q3,
                                                           rel=1e-12)


def test_divided_differences_recover_the_derivative():
    """`second_divided_difference` is q'''/2 and `third_divided_difference` q'''/6.

    The normalisations differ, and pairing an estimate with the wrong constant is
    exactly how the 3/4 optimism survived in two places, so both are pinned.
    Non-uniform grid again: the charge-based form is exact for a cubic on ANY grid
    and a formulation with a ratio-dependent bias would fail here.
    """
    q3 = 1.0e18
    q = lambda t: q3 * t ** 3 / 6.0
    g = lambda t: q3 * t ** 2 / 2.0           # g = dq/dt
    h1, h2, h3 = 1e-6, 6.1e-7, 3.3e-7
    t_n = 0.0

    dd2 = K.second_divided_difference(
        g(t_n), g(t_n - h1), g(t_n - h1 - h2), h1, h2)
    assert dd2 == pytest.approx(q3 / 2.0, rel=1e-6)

    dd3 = K.third_divided_difference(
        q(t_n), [q(t_n - h1), q(t_n - h1 - h2), q(t_n - h1 - h2 - h3)],
        h1, h2, h3)
    assert dd3 == pytest.approx(q3 / 6.0, rel=1e-6)


def test_both_backends_reach_the_same_definition():
    """The property the module exists for: one definition, two callers.

    Checked by identity rather than by value -- two copies that happen to agree
    today is the state this change replaced, and a value check would pass against
    it.  `integrator.py` re-exports the charge-based estimator under its old name;
    that must be the kernel's object, not a same-named twin.
    """
    from pycircuit.circuit import integrator as cpu
    assert cpu.third_divided_difference is K.third_divided_difference

    jax = pytest.importorskip('pycircuit.circuit.jaxtransient')
    ## Both backends' Gear-2 companion must BE the kernel call, so a future edit to
    ## one cannot drift from the other.
    assert jax.gear2_step(1.7, 2e-9, 1.3, 0.8, 1e-5, 7e-6) == \
        K.bdf2_companion(1.7, 2e-9, 1.3, 0.8, 1e-5, 7e-6)


def test_kernels_accept_arrays_elementwise():
    """Vectors, not just scalars -- the backends pass whole state vectors."""
    h1, h2 = 1e-5, 7e-6
    q_n = np.array([1.7, -0.3, 0.0])
    q_1 = np.array([1.3, -0.1, 0.0])
    q_2 = np.array([0.8, 0.2, 0.0])
    C = np.array([2e-9, 1e-12, 0.0])

    iq, geq = K.bdf2_companion(q_n, C, q_1, q_2, h1, h2)
    assert iq.shape == q_n.shape and geq.shape == C.shape
    for i in range(len(q_n)):
        s_iq, s_geq = K.bdf2_companion(q_n[i], C[i], q_1[i], q_2[i], h1, h2)
        assert iq[i] == pytest.approx(s_iq, rel=1e-15)
        assert geq[i] == pytest.approx(s_geq, rel=1e-15)
