"""Backend-agnostic scalar algebra for the integrators and their LTE estimates.

STAGE 9(a).  The formulas here were transcribed rather than shared, and every
defect stage 9 found was a repair that existed on one backend and not the other:

  * the Gear-2 error constant, fixed on the CPU in stage 4i and still 3/4
    optimistic on the JAX path until gate 9-1(a) measured it -- the same defect
    found three times in two copies;
  * the opening-step ramp, fixed on the CPU in stage 3 and missing from JAX until
    9(g), where it was pinning the error across four decades of ``reltol``;
  * the breakpoint scan, which 9(d) found wrong in two different ways in the two
    copies of the same loop, each hiding the other.

None of those was a new bug.  Each was a fix that did not cross the copy, which
is the argument for this module.

**These functions take no toolkit and import nothing.**  Every one is elementwise
arithmetic on charges, capacitances, companion currents and step sizes, so the
Python operators dispatch to whatever array type is passed -- ``numpy`` on the CPU
path, traced ``jax`` arrays inside ``jax.jit``.  That is not an accident of the
present implementation, it is the constraint: adding a ``numpy`` call here would
silently break tracing on the JAX backend, which is why gate 9a-1 asserts the
module imports neither backend.

Naming follows the transient literature and the rest of the package: ``h1`` is
the step about to be taken, ``h2`` the one before it, ``q`` a charge, ``g`` the
companion current ``dq/dt``.
"""


def bdf2_alphas(h1, h2):
    """The variable-step BDF-2 (Gear-2) coefficients.

    A second-order polynomial through the current point and the two before it,
    differentiated at the current point.  Fixed-step SPICE2 formulas are wrong the
    moment the step changes, which is the whole reason these are computed rather
    than tabulated.

    Returns ``(alpha0, alpha1, alpha2)`` such that
    ``dq/dt|_n ~ alpha0 q_n + alpha1 q_{n-1} + alpha2 q_{n-2}``.
    """
    return ((2.0 * h1 + h2) / (h1 * (h1 + h2)),
            -(h1 + h2) / (h1 * h2),
            h1 / (h2 * (h1 + h2)))


def bdf2_derivative(q_curr, q_prev1, q_prev2, h1, h2):
    """``dq/dt`` at the current point under Gear-2, without the conductance.

    Separate from :func:`bdf2_companion` because the LTE fallback needs to
    reconstruct ``g_n`` where no capacitance matrix is in scope; passing a dummy
    ``C`` to get one return value discarded is how a shared kernel starts growing
    arguments that mean nothing.
    """
    a0, a1, a2 = bdf2_alphas(h1, h2)
    return a0 * q_curr + a1 * q_prev1 + a2 * q_prev2


def bdf2_companion(q_curr, C_curr, q_prev1, q_prev2, h1, h2):
    """Gear-2 companion current and equivalent conductance.

    ``geq`` is the 'algorithmic conductance' the timestep introduces: the Jacobian
    contribution ``d(iq)/dx = alpha0 * C``.
    """
    a0, _a1, _a2 = bdf2_alphas(h1, h2)
    return (bdf2_derivative(q_curr, q_prev1, q_prev2, h1, h2), a0 * C_curr)


def euler_companion(q_curr, C_curr, q_prev, h):
    """Backward Euler companion current and equivalent conductance."""
    return (q_curr - q_prev) / h, C_curr / h


def trapezoidal_companion(q_curr, C_curr, q_prev, iq_prev, h):
    """Trapezoidal companion current and equivalent conductance.

    ``iq_n = 2 (q_n - q_{n-1})/h - iq_{n-1}`` -- note this recursion carries an
    undamped ``(-1)^n`` homogeneous mode, which is why the trapezoidal LTE
    estimator must NOT difference ``g``; see :func:`third_divided_difference` and
    stage 4g(b).
    """
    return 2.0 * (q_curr - q_prev) / h - iq_prev, 2.0 * C_curr / h


def second_divided_difference(g_n, g_nm1, g_nm2, h1, h2):
    """Second divided difference of the companion current: ``q'''/2``.

    ``g = dq/dt``, so differencing it twice gives ``g'' = q'''``, and the divided
    difference carries the conventional factor of ``1/2``.
    """
    return ((g_n - g_nm1) / h1 - (g_nm1 - g_nm2) / h2) / (h1 + h2)


def third_divided_difference(q_curr, q_last, h_curr, h_last, h_last2):
    """Third divided difference of the CHARGE: ``q'''/6``.

    STAGE 4g(b)/4i.  Taken from charges rather than from companion currents
    because the trapezoidal recursion above carries an undamped ``(-1)^n`` mode:
    an estimator that differences ``g`` measures that parasitic mode as well as
    the truncation error, and the mode does not shrink with ``h``.  Differencing
    the charge is mode-free and ratio-independent -- the mean-value form is exact
    for a cubic on any grid, uniform or not.

    ``q_last`` is the ring buffer, most recent first.
    """
    d1 = (q_curr - q_last[0]) / h_curr
    d2 = (q_last[0] - q_last[1]) / h_last
    d3 = (q_last[1] - q_last[2]) / h_last2

    dd_a = (d1 - d2) / (h_curr + h_last)
    dd_b = (d2 - d3) / (h_last + h_last2)

    return (dd_a - dd_b) / (h_curr + h_last + h_last2)


## The local truncation error of each method, given an estimate of the derivative
## it differences.  Kept next to the companion formulas deliberately: the pairing
## is what went wrong in 4i and again on the JAX path, where a GEAR2 residual
## taken from a different derivation reported (1/4) h^2 q''' against a true (1/3).

def euler_lte(g_n, g_nm1):
    """Backward Euler: ``-(1/2) h q''``, formed as a difference of companion
    currents.  A CURRENT, not a charge -- see ``Integrator.compute_lte``."""
    return -0.5 * (g_n - g_nm1)


def trapezoidal_lte(h1, q_triple_prime_over_6):
    """Trapezoidal: ``-(h^2/6) q'''``.

    Takes ``q'''/6`` -- the value :func:`third_divided_difference` returns -- so
    the caller cannot pair the formula with a differently-normalised estimate.
    """
    return -(h1 ** 2) * q_triple_prime_over_6


def gear2_lte(h1, h2, q_triple_prime_over_6):
    """Gear-2: ``-(1/6) h1 (h1 + h2) q'''``.

    On a uniform grid this is ``-(1/3) h^2 q'''``.  YWR's Table I GEAR2 residual
    gives ``(1/4)`` instead, i.e. 3/4 of the truth; that formula was used on the
    CPU until stage 4i and on the JAX path until gate 9-1(a), so the constant is
    derived here once rather than transcribed twice.
    """
    return -h1 * (h1 + h2) * q_triple_prime_over_6
