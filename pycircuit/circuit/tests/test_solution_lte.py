"""STAGE 12B -- Fang eq (6), the solution-space truncation error.

    eps_m = | v_i(t_m) - v_{i,extrapolated} |

Stage 12 was built for a long time against a charge divided-difference LTE, which
is a different estimator with the same order.  These tests pin the properties
that made the substitution matter, so a future edit cannot quietly swap it back:

  * the extrapolation is EXACT for polynomials up to its degree -- that is what
    makes the deviation an error estimate rather than a difference of two
    approximations;
  * the deviation converges as h^(k+1), so it measures truncation and not noise;
  * and, the property the whole coupled method rests on, it does NOT blow up as
    h shrinks, which the charge-based estimator does.

The last one is tested by comparison against the existing kernel, because a
statement about conditioning is only meaningful relative to something.
"""
import numpy as np
import pytest

from pycircuit.circuit import _lte_kernels as K


def test_extrapolation_times_run_backwards_from_the_last_point():
    """Origin at the most recent accepted point, past negative.

    Keeps the target `+h_curr` a small positive number instead of a difference of
    large absolute times, which is where an adaptive run with `t` in the seconds
    and `h` in the picoseconds loses its significant digits.
    """
    assert K.extrapolation_times([2.0, 3.0, 5.0]) == [0.0, -2.0, -5.0, -10.0]
    assert K.extrapolation_times([]) == [0.0]


@pytest.mark.parametrize('degree', [0, 1, 2, 3])
def test_extrapolation_is_exact_for_polynomials_up_to_its_degree(degree):
    """`degree+1` points reproduce a degree-`degree` polynomial exactly.

    On a deliberately NON-uniform grid: equally spaced points would let a
    fixed-step formula pass, and an adaptive run never produces one.
    """
    coef = [0.3, -1.7, 2.1, -0.9][:degree + 1]
    poly = lambda t: sum(c * t ** i for i, c in enumerate(coef))

    h_hist = [1e-5, 3.7e-6, 8.1e-6]
    times = K.extrapolation_times(h_hist)[:degree + 1]
    v_hist = [poly(t) for t in times]

    h_curr = 2.3e-6
    got = K.extrapolate(v_hist, h_hist, h_curr)
    assert got == pytest.approx(poly(h_curr), rel=1e-9, abs=1e-18)


def test_deviation_vanishes_when_the_solution_is_the_fitted_polynomial():
    """If the circuit followed the polynomial exactly, the LTE is zero."""
    h_hist = [1e-5, 4e-6]
    poly = lambda t: 2.0 - 3.0 * t + 7.0 * t ** 2
    times = K.extrapolation_times(h_hist)[:3]
    v_hist = [poly(t) for t in times]
    h_curr = 6e-6
    lte = K.solution_lte(poly(h_curr), v_hist, h_hist, h_curr)
    assert abs(lte) < 1e-15


def test_deviation_converges_at_the_order_of_the_extrapolation():
    """Degree-k extrapolation must leave an O(h^(k+1)) deviation.

    Measured as an observed exponent from two step sizes a decade apart rather
    than asserted against a constant, so the test states the ORDER rather than a
    fitted number.  `exp(t)` is used because every derivative is nonzero, so no
    term of the error expansion vanishes by accident.
    """
    for degree, expect in ((1, 2.0), (2, 3.0)):
        errs = []
        for h in (1e-3, 1e-4):
            ## Uniform past grid here: the quantity under test is the exponent,
            ## and a varying grid would change the constant between the two runs.
            h_hist = [h] * degree
            times = K.extrapolation_times(h_hist)[:degree + 1]
            v_hist = [np.exp(t) for t in times]
            errs.append(abs(K.solution_lte(np.exp(h), v_hist, h_hist, h)))
        observed = np.log10(errs[0] / errs[1])
        assert observed == pytest.approx(expect, abs=0.15), \
            'degree %d: observed order %.3f, expected %.1f' % (degree, observed, expect)


def test_it_is_elementwise_over_a_state_vector():
    """Whole solution vectors, and each unknown independent of the others."""
    h_hist = [1e-5, 4e-6]
    times = K.extrapolation_times(h_hist)[:3]
    fns = (lambda t: 1.0 + t, lambda t: -2.0 * t ** 2, lambda t: 0.0 * t + 3.0)
    v_hist = [np.array([f(t) for f in fns]) for t in times]
    h_curr = 6e-6

    got = K.extrapolate(v_hist, h_hist, h_curr)
    assert got.shape == (3,)
    for i, f in enumerate(fns):
        assert got[i] == pytest.approx(f(h_curr), rel=1e-9, abs=1e-18)


def test_it_does_not_blow_up_as_the_step_shrinks():
    """THE PROPERTY THE COUPLED METHOD RESTS ON.

    Fang solves `f_lte(v, h) = 0` for `h`, so the estimator has to stay
    well-behaved as Newton moves `h` around.  The charge-based estimator does
    not: `third_divided_difference` divides by `h_curr` twice, so with a fixed
    perturbation in the stored charges its output grows without bound as
    `h_curr -> 0`, and a Newton iteration chasing it walks the step size down
    until it underflows -- which is exactly what stage 12B-0 measured.

    Both estimators are given the same rounding-level perturbation and swept over
    four decades of `h_curr`.  The solution-space form must stay bounded; the
    charge form is asserted to diverge, so this test also documents WHY the
    estimator was changed and fails if that ever stops being true.
    """
    noise = 1e-16
    hs = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    solution_vals, charge_vals = [], []
    for h in hs:
        h_hist = [1e-6, 1e-6, 1e-6]
        times = K.extrapolation_times(h_hist)
        ## A perfectly smooth history, plus one unit of rounding on the newest
        ## point -- the situation at a converged solve.
        v_hist = [np.cos(1e5 * t) for t in times[:3]]
        v_curr = np.cos(1e5 * h) + noise
        solution_vals.append(abs(K.solution_lte(v_curr, v_hist, h_hist, h)))

        q_last = [np.cos(1e5 * t) for t in times[1:4]]
        q_curr = np.cos(1e5 * h) + noise
        charge_vals.append(abs(K.third_divided_difference(
            q_curr, q_last, h, 1e-6, 1e-6)))

    ## What matters is the DIRECTION, not the magnitude: over this sweep the
    ## history spacing is held fixed while `h_curr` shrinks, so neither value is
    ## a clean convergence rate -- but one must fall and the other must rise.
    ## Asserting an absolute floor here was wrong: the solution-space value at
    ## h=1e-6 is genuine truncation error (~5e-5 on this signal), not noise.
    sol_ratio = solution_vals[-1] / solution_vals[0]
    chg_ratio = charge_vals[-1] / charge_vals[0]

    assert sol_ratio < 1e-3, \
        'solution-space LTE must fall as h shrinks, ratio was %g' % sol_ratio

    ## The charge-based one goes the other way, and by orders of magnitude. That
    ## is the disqualifying property: Newton solving f_lte = 0 for h then walks
    ## the step down forever, which is what stage 12B-0 measured.
    assert chg_ratio > 1e3, \
        'charge divided difference no longer diverges as h->0 (ratio %g); ' \
        'if this is now false, revisit doc/fang_dac2013_math.md sec. 3' % chg_ratio

    ## No third assertion on the gap between the two ratios: the obvious one is a
    ## threshold with nothing behind it, and a bound picked after seeing the
    ## measurement is fitted rather than derived.  The two directional assertions
    ## above are the property; the sweep is printed by the 12B benchmark for
    ## anyone who wants the magnitudes.
