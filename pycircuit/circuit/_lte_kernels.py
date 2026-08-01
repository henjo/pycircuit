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


## ---------------------------------------------------------------------------
## STAGE 12B -- the SOLUTION-SPACE truncation error, Fang DAC 2013 eq (6).
##
##     eps_m = | v_i(t_m) - v_{i,extrapolated} |
##
## "a typical industrial circuit simulator calculates the difference between the
## computed solution and the polynomial extrapolation from previous time steps,
## and considers the maximum value a good estimate for the local truncation
## error [...] the corresponding node is often referred to as controlling LTE
## node, which may vary from time point to time point."
##
## This is a DIFFERENT ESTIMATOR from the charge-based ones above, not a
## reformulation of them, and the difference is the point.  The charge divided
## differences carry repeated `1/h` factors, so as the step shrinks they amplify
## rounding in `q` and the estimate GROWS -- usable as a test, useless as an
## equation to solve for `h`, because Newton then runs away from the root.  The
## form below differences two solution values and has no such amplification,
## which is what makes Fang's coupled system solvable at all.
##
## See `doc/fang_dac2013_math.md` sec. 3.


def extrapolation_times(h_hist):
    """Times of the stored past points, relative to the most recent one.

    ``h_hist[j]`` is the step that ENDED at past point ``j``, so ``h_hist[0]`` is
    the step from point 1 to point 0.  Returns ``[0, -h_hist[0],
    -(h_hist[0]+h_hist[1]), ...]`` -- the most recent accepted point is the
    origin and the past runs backwards, which keeps the target time ``+h_curr``
    a small positive number rather than a difference of large absolutes.
    """
    times = [0.0]
    acc = 0.0
    for h in h_hist:
        acc = acc - h
        times.append(acc)
    return times


def extrapolate(v_hist, h_hist, h_curr):
    """Polynomial extrapolation of the solution to ``t_{m-1} + h_curr``.

    ``v_hist`` is the accepted solution history, most recent first; ``k+1``
    entries give a degree-``k`` polynomial.  Newton's divided-difference form is
    used rather than a fitted Vandermonde system because the past points are
    NOT equally spaced -- an adaptive run never produces a uniform grid, and a
    fixed-step formula is wrong the moment the step changes.

    Elementwise, so ``v_hist`` may be scalars or whole state vectors.
    """
    n = len(v_hist)
    if n == 0:
        raise ValueError('extrapolate needs at least one past point')

    times = extrapolation_times(h_hist)[:n]
    if len(times) < n:
        raise ValueError('need %d step sizes for %d points, got %d'
                         % (n - 1, n, len(h_hist)))

    ## Divided-difference table, built on a copy so the caller's history is not
    ## consumed.  `dd[i]` descends to the i-th order difference in place.
    dd = list(v_hist)
    coeffs = [dd[0]]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / (times[i] - times[i - j])
        coeffs.append(dd[j])

    ## Horner-like accumulation of c0 + c1 (t-t0) + c2 (t-t0)(t-t1) + ...
    out = coeffs[0]
    basis = 1.0
    for j in range(1, n):
        basis = basis * (h_curr - times[j - 1])
        out = out + coeffs[j] * basis
    return out


def extrapolation_error_weight(h_curr, h_hist):
    """The node polynomial that scales the extrapolation deviation.

    For a degree-``k`` polynomial through the accepted points, the interpolation
    error at ``t = h_curr`` is

        E(h) = (v^(k+1)(xi) / (k+1)!) * h (h+h1) (h+h1+h2) ... (h+h1+..+hk)

    so this returns the product, which is everything in the deviation that
    depends on the step size.

    **It is NOT proportional to h^(k+1) except in one limit.**  When
    ``h >> h1`` the product behaves like ``h^(k+1)``; when ``h << h1`` every
    factor but the first is roughly constant and the deviation is LINEAR in
    ``h``.  A controller that assumes the ``h^(k+1)`` law everywhere therefore
    mis-predicts by orders of magnitude on the shrinking side -- measured on
    `rc-vsin` at reltol 1e-6 as 2049 rejections against 2143 accepted steps,
    oscillating between a step far under tolerance and one far over.

    This is also Fang's ``d = df_lte/dh_m`` up to the constant: differentiating
    the product is closed form, which is what lets sec. 3.2 claim ``p, q^T`` and
    ``d`` cost next to nothing.
    """
    weight = h_curr
    acc = 0.0
    for h in h_hist:
        acc = acc + h
        weight = weight * (h_curr + acc)
    return weight


def step_for_error_ratio(h_curr, h_hist, ratio, lo_factor, hi_factor, iters=40):
    """Step size whose extrapolation weight is ``ratio`` times the current one.

    Inverts :func:`extrapolation_error_weight`, which is strictly increasing in
    ``h_curr`` for positive steps, so a bisection on ``[lo_factor*h_curr,
    hi_factor*h_curr]`` is monotone and always terminates.  The bracket is the
    controller's own growth/shrink clamp, so the search never proposes a step it
    would have rejected anyway, and no separate clamp is needed afterwards.

    A closed-form cubic root would be marginally cheaper for degree 2 and would
    not generalise to another order; this costs a few dozen scalar operations
    once per step, against a Newton solve of the whole circuit.
    """
    target = extrapolation_error_weight(h_curr, h_hist) * ratio

    lo, hi = h_curr * lo_factor, h_curr * hi_factor
    if extrapolation_error_weight(lo, h_hist) >= target:
        return lo
    if extrapolation_error_weight(hi, h_hist) <= target:
        return hi

    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if extrapolation_error_weight(mid, h_hist) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solution_lte(v_curr, v_hist, h_hist, h_curr):
    """Fang eq (6) before the max: the per-unknown extrapolation deviation.

    Returns ``v_curr - P(t_m)`` elementwise.  Taking the maximum, and deciding
    what to normalise by first, is the caller's business -- eq (6) maximises the
    raw difference, but a vector mixing volts and amps has to be normalised by a
    per-unknown tolerance before a maximum means anything, and the paper does not
    specify ``tau_m``.  That choice is recorded in `doc/fang_dac2013_math.md`
    sec. 6 as ours rather than Fang's.
    """
    return v_curr - extrapolate(v_hist, h_hist, h_curr)


def gear2_lte(h1, h2, q_triple_prime_over_6):
    """Gear-2: ``-(1/6) h1 (h1 + h2) q'''``.

    On a uniform grid this is ``-(1/3) h^2 q'''``.  YWR's Table I GEAR2 residual
    gives ``(1/4)`` instead, i.e. 3/4 of the truth; that formula was used on the
    CPU until stage 4i and on the JAX path until gate 9-1(a), so the constant is
    derived here once rather than transcribed twice.
    """
    return -h1 * (h1 + h2) * q_triple_prime_over_6
