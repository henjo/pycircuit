"""STAGE 12B -- ``d(iq)/dh``, the integrator half of Fang's ``p = df_ckt/dh``.

Once the step size is an unknown (DAC 2013 eq 11) the circuit residual has to be
differentiated with respect to it, and it depends on ``h`` through two routes:
the sources at ``t_{m-1}+h`` (``TimeFunction.dfdt``, tested separately) and the
companion current. The paper writes the integration coefficients as explicit
functions of ``h_m`` in eq (4) for exactly this reason.

Every derivative here is checked against a central difference of the companion
current it differentiates, on a NON-UNIFORM history -- a uniform grid lets a
fixed-step formula pass, and an adaptive run never produces one.
"""
import numpy as np
import pytest

from pycircuit.circuit import _lte_kernels as K
from pycircuit.circuit.integrator import (EulerIntegrator, TrapezoidalIntegrator,
                                          Gear2Integrator)
from pycircuit.circuit.toolkit import numeric


H1, H2 = 1.3e-6, 8.0e-7
Q_CURR = np.array([1.7e-9, -3.1e-10, 5.0e-9])
Q_1 = np.array([1.5e-9, -2.6e-10, 4.2e-9])
Q_2 = np.array([1.1e-9, -1.9e-10, 3.7e-9])
C_CURR = np.array([2e-12, 1e-13, 8e-13])
IQ_1 = np.array([1.4e-4, -5.0e-5, 3.3e-4])


def _companion(integ, h):
    """The companion CURRENT only -- the conductance is a Jacobian term."""
    iq, _geq = integ.compute_derivatives(
        q_curr=Q_CURR, C_curr=C_CURR, h_curr=h,
        q_last=[Q_1, Q_2], iq_last=[IQ_1, IQ_1], h_last=H2,
        is_first_step=False, toolkit=numeric)
    return iq


CASES = [
    ('euler', EulerIntegrator()),
    ('trapezoidal', TrapezoidalIntegrator()),
    ('gear2', Gear2Integrator()),
]


@pytest.mark.parametrize('name,integ', CASES, ids=[c[0] for c in CASES])
def test_companion_dh_matches_a_central_difference(name, integ):
    """The derivative against the function it differentiates."""
    analytic = integ.companion_dh(Q_CURR, [Q_1, Q_2], H1, H2)
    eps = 1e-6 * H1
    fd = (_companion(integ, H1 + eps) - _companion(integ, H1 - eps)) / (2 * eps)
    assert np.allclose(analytic, fd, rtol=1e-5), \
        '%s: analytic %r, finite difference %r' % (name, analytic, fd)
    ## Not vacuously zero: a gradient test passing on 0 == 0 checks nothing.
    assert np.max(np.abs(analytic)) > 0.0


def test_bdf2_alpha_derivatives_match_a_difference_of_the_alphas():
    """`bdf2_alphas_dh` against `bdf2_alphas`, coefficient by coefficient.

    Separate from the companion check because a sign error in one alpha can be
    masked by the charges it multiplies -- with three charges of the same sign,
    two compensating errors give a plausible total.
    """
    eps = 1e-7 * H1
    analytic = K.bdf2_alphas_dh(H1, H2)
    ap = K.bdf2_alphas(H1 + eps, H2)
    am = K.bdf2_alphas(H1 - eps, H2)
    for i in range(3):
        fd = (ap[i] - am[i]) / (2 * eps)
        assert analytic[i] == pytest.approx(fd, rel=1e-5), \
            'alpha%d: analytic %g, finite difference %g' % (i, analytic[i], fd)


def test_only_the_current_step_is_a_variable():
    """`h2` is a past step and is held fixed.

    Differentiating with respect to it as well would be differentiating a
    constant, and would put a term in `p` for a quantity the coupled system does
    not solve for.
    """
    d_at_h2 = K.bdf2_alphas_dh(H1, H2)
    d_at_other_h2 = K.bdf2_alphas_dh(H1, 2.0 * H2)
    ## Sanity: the derivative depends on h2 as a parameter, but nothing here
    ## differentiates with respect to it -- the two differ, which is expected,
    ## and neither is a derivative in h2.
    assert d_at_h2 != d_at_other_h2


def test_an_integrator_without_the_derivative_says_so():
    """A silent zero would look like a converged solve of the wrong problem."""
    from pycircuit.circuit.integrator import Integrator

    class Bare(Integrator):
        ORDER = 1
        def check_order_drop(self, h_curr, h_last, is_first_step):
            return self
        def compute_derivatives(self, *a, **k):
            return None, None
        def compute_lte(self, *a, **k):
            return None, 1.0
        def get_required_history(self):
            return 1

    with pytest.raises(NotImplementedError) as exc:
        Bare().companion_dh(Q_CURR, [Q_1, Q_2], H1, H2)
    assert 'coupled' in str(exc.value)
