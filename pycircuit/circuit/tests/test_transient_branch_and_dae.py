"""Shooting tests: transient branch and dae.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
from pycircuit.circuit.tests._shooting_fixtures import (_comparator_relaxation_oscillator,
    _q20_rlc,
    _rc_ladder)


def test_the_alternating_mode_is_l_stability_not_the_iq_recursion():
    """Why trapezoidal shooting dies at even step counts -- the real cause.

    The class docstring blamed trapezoidal's `iq` RECURSION
    (`iq_n = ... - iq_{n-1}`, homogeneous mode `(-1)^n`).  That is the wrong
    cause: trapezoidal is A-stable but NOT L-stable, so it maps the null
    space of the singular MNA `C` by exactly `-1` per step, and an X-ONLY
    formulation carrying no `iq` variable at all is singular in the same
    way.  For `C x' + G x + u = 0`,

        A_trap  = (C/h + G/2)^-1 (C/h - G/2)   ->  -I  on null(C)
        A_euler = (C/h + G)^-1 (C/h)           ->   0  on null(C)

    so the count of `-1` modes is exactly `m - rank(C)`, one per ALGEBRAIC
    variable, on every MNA circuit.  That is a prediction with a number in
    it, which is why it is worth a test: it says the opening step must be
    L-STABLE and nothing more specific, so any L-stable opener rescues
    exactly those modes and Euler is not privileged.
    """
    circuit.default_toolkit = circuit.numeric
    h = 1e-3 / 200
    for name, cir in (('Q=20 RLC', _q20_rlc()), ('RC', _rc_ladder(1))):
        n = cir.n
        Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)
        Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)
        iref = cir.get_node_index(gnd)
        keep = [i for i in range(n) if i != iref]
        Cm = Cm[np.ix_(keep, keep)]
        Gm = Gm[np.ix_(keep, keep)]
        m = len(keep)
        rank_c = np.linalg.matrix_rank(Cm)
        algebraic = m - rank_c
        assert algebraic > 0, '%s has no algebraic variables to test' % name

        e_trap = np.linalg.eigvals(np.linalg.solve(Cm / h + Gm / 2.0,
                                                   Cm / h - Gm / 2.0))
        e_eul = np.linalg.eigvals(np.linalg.solve(Cm / h + Gm, Cm / h))
        n_minus1 = int(np.sum(np.abs(e_trap + 1.0) < 1e-9))
        n_zero = int(np.sum(np.abs(e_eul) < 1e-9))
        assert n_minus1 == algebraic, \
            '%s: trapezoidal has %d modes at -1, m - rank(C) = %d' \
            % (name, n_minus1, algebraic)
        assert n_zero == algebraic, \
            '%s: Euler has %d modes at 0, m - rank(C) = %d' \
            % (name, n_zero, algebraic)


def test_no_limiter_in_the_tree_has_a_charge_that_reads_its_limiting_state():
    """`_C_at` carries NO limiting sync, and this is the measurement that says
    it may not need one.

    `_G_at` must sync: a junction's `i`/`G` are read at the device's stored
    `_vlim`, which is what the shared-`_vlim` monodromy defect was about.  CHARGE
    is a different question, and the answer across this tree is that no device's
    `C`/`q` reads that state:

      * `elements.Diode` is the ONLY stateful limiter (it keeps `_vlim`);
      * `Semiconductor` (BJT/JFET/ZenerDiode/Varactor) limits STATE-FREE by
        construction, as its own docstring says;
      * `compact.PspMosLongChannel` likewise returns a limited copy;
      * the hdl devices keep no `_vlim` (it is a code-generation local).

    So the only device that COULD show the effect is the plain `Diode`, and it
    does not.

    ⚠ THE CONTROL IS THE POINT OF THIS TEST.  `dC = 0` on its own is worthless:
    it is equally what you get if `limit()` did nothing at all, which is exactly
    what happened when this was first probed with two hdl devices that do not
    respond to `cir.limit` (`dG = di = 0`, a vacuous pass).  So the assertion
    below REQUIRES the conductance to move -- proving the limiting really was
    live -- before it accepts that the charge did not.

    If a stateful limiter whose charge DOES read its state is ever added, this
    test fails and `_C_at` needs its sync back.
    """
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric
    ep = circuit.defaultepar.copy()

    def build():
        c = SubCircuit()
        c['vs'] = VS(1, gnd, v=0.0)
        c['R'] = R(1, 2, r=1e3)
        c['D'] = Diode(2, gnd)
        return c

    c1 = build()
    k = c1.get_node_index(2)
    x = np.zeros(c1.n)
    x[k] = 0.75                      # the junction well into conduction
    C0 = np.asarray(c1.C(x, ep), dtype=float).copy()
    q0 = np.asarray(c1.q(x, ep), dtype=float).copy()
    G0 = np.asarray(c1.G(x, ep), dtype=float).copy()
    i0 = np.asarray(c1.i(x, ep), dtype=float).copy()

    ## the same point, but with the device's limiting state left far away
    c2 = build()
    xs = np.zeros(c2.n)
    xs[c2.get_node_index(2)] = 0.05
    c2.limit(xs, xs, ep)
    C1 = np.asarray(c2.C(x, ep), dtype=float)
    q1 = np.asarray(c2.q(x, ep), dtype=float)
    G1 = np.asarray(c2.G(x, ep), dtype=float)
    i1 = np.asarray(c2.i(x, ep), dtype=float)

    ## CONTROL FIRST: the limiting must actually have taken effect, or a null
    ## result below means nothing at all.
    dG = float(np.max(np.abs(G1 - G0)))
    di = float(np.max(np.abs(i1 - i0)))
    assert dG > 1.0 and di > 1e-3, \
        'the control failed: moving the stored limiting state changed G by ' \
        '%.2e and i by %.2e, so `limit()` did nothing here and the charge ' \
        'result below would be vacuous' % (dG, di)

    dC = float(np.max(np.abs(C1 - C0)))
    dq = float(np.max(np.abs(q1 - q0)))
    assert dC == 0.0 and dq == 0.0, \
        'a limiter whose CHARGE reads its stored state now exists (dC=%.2e, ' \
        'dq=%.2e) -- `_C_at` needs its `_sync_limit_at` back' % (dC, dq)

def test_fsolve_stops_a_stalled_iteration_when_asked():
    """`fsolve(stall_window=W)` (2026-09-26): stop when the best `||F||` of
    the last `W` iterations is not half the best before them, with `ier = 2`
    and `infodict['stalled']`.  Off by default.  `x^2 + 1` has no real root:
    Newton wanders, and stops after W + 1 evaluations instead of `maxiter`."""
    import numpy as _np
    import pycircuit.circuit.analysis as _an
    from pycircuit.circuit import numeric as _tk
    n = [0]

    def f(x):
        n[0] += 1
        return _np.array([x[0] ** 2 + 1.0]), _np.array([[2.0 * x[0] + 1e-12]])
    _x, info, ier, mesg = _an.fsolve(f, _np.array([0.7]), maxiter=200,
                                     toolkit=_tk, full_output=True,
                                     stall_window=5)
    assert ier == 2 and info['stalled'] and 'Stalled' in mesg, (ier, info, mesg)
    assert n[0] < 30, n[0]
    n[0] = 0
    _x, info, ier, _m = _an.fsolve(f, _np.array([0.7]), maxiter=200,
                                   toolkit=_tk, full_output=True)
    assert ier == 2 and not info['stalled'] and n[0] == 200, (info, n[0])


def test_a_line_search_trial_that_cannot_be_evaluated_is_halved_not_fatal():
    """`fsolve(line_search=True)` (2026-09-26): a trial whose evaluation
    raises `NoConvergenceError` -- a shooting residual stepping a transient
    from a state far off the orbit -- counts as uphill and is halved.  It
    used to abort the whole solve from a trial the halving would have pulled
    back (gear's free-period solve on a smooth 3:1 grid died that way at 200
    and 240 points).  Newton on arctan from 1.5 overshoots to -1.69; the
    residual refuses |x| > 1.6: the solve now converges.  When NO trial in
    the budget evaluates, the error still propagates."""
    import numpy as _np
    import pycircuit.circuit.analysis as _an
    from pycircuit.circuit import numeric as _tk

    def f(x):
        v = float(x[0])
        if abs(v) > 1.6:
            raise _an.NoConvergenceError('the step diverged')
        return (_np.array([_np.arctan(v)]),
                _np.array([[1.0 / (1.0 + v * v)]]))
    x, _i, ier, _m = _an.fsolve(f, _np.array([1.5]), maxiter=40,
                                toolkit=_tk, full_output=True,
                                line_search=True)
    assert ier == 1 and abs(float(x[0])) < 1e-8, (ier, x)

    def g(x):
        if abs(float(x[0]) - 1.5) > 1e-9:
            raise _an.NoConvergenceError('every trial diverges')
        return (_np.array([_np.arctan(float(x[0]))]),
                _np.array([[1.0 / (1.0 + float(x[0]) ** 2)]]))
    with pytest.raises(_an.NoConvergenceError, match='every trial'):
        _an.fsolve(g, _np.array([1.5]), maxiter=40, toolkit=_tk,
                   full_output=True, line_search=True)



def test_the_sources_off_numerical_floor_drops_all_the_way_to_the_arithmetic():
    """⚠⚠ THE SOURCES-OFF NOISE FLOOR, measured on THIS tree's construction
    rather than on the published one.

    Biggio/Bizzarri/Brambilla/Storace 2013 measure the floor by FFT-ing the
    jitter of a time-domain noise run with the sources off.  That does NOT
    transfer here: with the sources off, `c` is identically zero and the
    closed-form PSD is exactly zero -- a structural identity that cannot fail,
    already recorded in the roadmap as a rejected gate.  What DOES transfer is
    the underlying object: threshold crossings of a simulated waveform, whose
    spacing carries the integrator's period error with or without a source.

    So the floor is measured the way a floor is measured on a bench -- inject a
    KNOWN deterministic perturbation of amplitude `a`, require the estimator to
    track it with slope 1, then turn it off and read what is left.  The sloped
    region is what makes this impossible to pass zero-versus-zero.

    Asserted, all measured (van der Pol Q = 15.9, 120 points per period, 100
    cycles with 80 discarded):

    1. the estimator is LINEAR in the injected perturbation, slope ~1 -- with
       no sloped region a plateau means nothing;
    2. radau's floor is far below gear-2's on the SAME grid.  This is
       docs-46's falsifiable prediction from the paper's order argument, and
       it holds by 2891x here (5.3e-10 against 2.3e-13);
    3. radau's floor is at the FLOATING-POINT limit of the time variable --
       measured 13x `ulp(t)/T0` -- so the order lever does not merely lower
       this floor, it reaches the bottom of the arithmetic, where no further
       order or grid refinement can help;
    4. the floor is DETERMINISTIC, not noise: a repeat run reproduces the
       period sequence to the bit.  A deterministic period error is a
       FREQUENCY SHIFT and makes no jitter at all; what makes the variation is
       the crossing landing somewhere different inside a step each cycle.

    ⚠⚠ FIVE INSTRUMENT FAILURES PRECEDED THESE NUMBERS, four of them caught by
    one tell -- TWO DIFFERENT INTEGRATORS AGREEING TO FIVE SIGNIFICANT FIGURES
    on a quantity that is supposed to BE their own error.  They are recorded in
    `benchmarks/noise_floor_sources_off.py`; the one that survived fixing the
    other four was the FIXTURE's own settling (tau = 31.8 cycles, and 5 were
    being discarded).  The full order-and-grid sweep lives in that benchmark
    because it costs minutes; this pins the instrument that produced it and the
    one comparison the result rests on.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.integrator import (Gear2Integrator,
                                              RadauIIA3Integrator)
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from noise_floor_sources_off import periods

    NPTS, NCYC, DROP = 120, 100, 80

    def jitter(cls, a):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            Tk, T0, _ = periods(cls, NPTS, a=a, ncyc=NCYC, drop=DROP)
        return float(np.std(Tk)) / T0, Tk, T0

    ## (1) the instrument is alive and linear
    j3, _, _ = jitter(RadauIIA3Integrator, 1e-3)
    j4, _, _ = jitter(RadauIIA3Integrator, 1e-4)
    slope = np.log(j3 / j4) / np.log(10.0)
    assert abs(slope - 1.0) < 0.05, (slope, j3, j4)

    ## (2) the order lever, on the same grid.  Measured 2891x; gated at 100x.
    floor_gear, _, _ = jitter(Gear2Integrator, 0.0)
    floor_radau, Tk1, T0 = jitter(RadauIIA3Integrator, 0.0)
    assert floor_radau < floor_gear / 100.0, (floor_gear, floor_radau)

    ## (3) and radau is at the arithmetic, not at its discretisation
    ulp_rel = float(np.spacing(NCYC * T0)) / T0
    assert floor_radau < 100.0 * ulp_rel, (floor_radau, ulp_rel)
    ## gear is NOT -- or (2) would be comparing two representation limits
    assert floor_gear > 1000.0 * ulp_rel, (floor_gear, ulp_rel)

    ## (4) deterministic, to the bit
    _, Tk2, _ = jitter(RadauIIA3Integrator, 0.0)
    assert np.array_equal(Tk1, Tk2)


def test_violating_the_im_D_hypothesis_costs_the_high_order_methods_their_order():
    """⚠⚠ THE `im D(t)` HYPOTHESIS BITES, AND IT COSTS RADAU EIGHT ORDERS.

    Lamour, März & Tischendorf (2013) put one hypothesis under IRK(DAE)
    convergence (Thm 5.7), GLM convergence at stage order (Thm 5.9 -- what
    `NordsieckGLMIntegrator` rests on) and contractivity transfer (Thm 6.9):
    `im D(t)` time-invariant.  It is one term in one equation -- the IERODE's
    field is `u' = R'(t)u + D(t)omega(u,t)` and the hypothesis exists to kill
    `R'(t)u`.  For charge-oriented MNA it means `im C(x)` constant along the
    orbit.

    ⚠ EVERY OTHER FIXTURE IN THIS TREE SATISFIES IT VACUOUSLY.  Measured
    2026-09-10: on the index-2 C-V loop, the state-free exponential and the van
    der Pol, `C(x)` is LITERALLY CONSTANT (`max|C(x1) - C(x2)| == 0`), because
    every reactance in them is linear.  And the bias is not particular to
    them -- a smoothly varying `C(v) > 0` is rank-1 throughout, so ordinary
    circuit fixtures cannot exercise the condition either.

    THE INSTRUMENT IS NOT AN ORDER SWEEP.  Example 3.34 / Thm 3.53: what fails
    at a regularity boundary is UNIQUENESS -- two solutions through a critical
    point -- not accuracy.  So this integrates through the crossing several
    times, differing ONLY in where the crossing lands inside a step (the grid
    is shifted; both endpoints stay pinned), and asks how that spread behaves
    under refinement.  ⚠ Shifting `tend` instead is an artefact and was this
    measurement's sixth instrument failure: the runs then end at DIFFERENT
    TIMES and the spread is `|dV/dt|*h*doffset`, which reads the same for
    every method and shrinks at `O(h)` for all of them.

    MEASURED, spread at 800 points per period:

    ======  ==========  ==============  ============
    method  LINEAR C    NONLINEAR C>0   rank C DROPS
    ======  ==========  ==============  ============
    gear    3.00e-06    2.93e-06        1.91e-05
    trbdf2  5.10e-09    7.55e-09        8.63e-06
    glm3    2.26e-12    3.65e-12        1.46e-05
    radau   4.44e-16    6.66e-16        9.75e-06
    ======  ==========  ==============  ============

    On both controls the spread shrinks at the method's own LOCAL order and
    the methods separate by ten orders; where rank C drops they collapse onto
    one magnitude and one slow rate.

    ⚠⚠ THE NONLINEAR CONTROL IS WHAT MAKES THAT MEAN ANYTHING -- without it
    the comparison confounds the rank change with the nonlinearity, since the
    violating fixture's `C` is nonlinear and the linear control's is not.  A
    `C(v) > 0` varying 2x across the orbit behaves IDENTICALLY to a constant
    one, so the loss is the RANK CHANGE.  ⚠ And `q = c0 V^3/3` is a
    polynomial, so it is not a smoothness failure of the model either.

    ⚠⚠ AND THE ATTRIBUTION IS NARROWER THAN IT LOOKS.  The numbers above are
    what they are, but this is NOT a DAE or index effect: the same collapse
    reproduces in a ONE-NODE SCALAR model with no DAE structure (docs-46,
    relayed), and its exponent is set by the ORDER OF THE ZERO OF `C`, not by
    the method's order.  MEASURED here on this fixture, fitting `spread ~ h^p`
    over N = 200/400/800/1600: with `C ~ |V|` (a first-order zero) TR-BDF2
    reads 1.496 and GLM3 1.527 against a predicted 1.5; with `C ~ V^2` they
    read 1.413 and 1.411 against a predicted 1.333.  radau is noisier (1.611
    and 1.709) and moves the wrong way with `k`, which is not explained.  So
    `im D(t)` is the right DESCRIPTION of when this happens -- rank C drops --
    without being the mechanism.

    ⚠ NOT SHOWN: the spread still shrinks, so this is an ORDER COLLAPSE and
    not the uniqueness failure the theory points at -- and per Lamour §2.9
    that is the PUBLISHED SIGNATURE OF A HARMLESS CRITICAL POINT (one "which
    disappears in smoother settings"), so it is the expected outcome here
    rather than a falsifier that missed.  ⚠⚠ This probe COULD NOT have seen a
    uniqueness failure anyway: on a genuine branch point the non-uniqueness
    becomes multiplicity of roots of the STEP equation, the Newton picks one
    silently, every grid picks the same one, and the spread over grid offset
    is exactly zero (relayed, measured by docs-46).  The knob that exposes it
    is the SOLVER'S INITIAL GUESS.  A real non-uniqueness fixture needs the
    rank drop AND a REPELLING equilibrium at it -- a negative conductance,
    which is what an oscillator's active device supplies.  NOT BUILT.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from im_d_falsifier import (rank_changing, constant_rank,
                                nonlinear_constant_rank, rank_probe, endpoint)

    ## (0) the instrument check: the fixtures must actually differ in rank, or
    ## nothing below is about the hypothesis at all
    ranks_bad = [r for _v, _m, r in rank_probe(rank_changing)]
    ranks_ok = [r for _v, _m, r in rank_probe(nonlinear_constant_rank)]
    assert min(ranks_bad) == 0 and max(ranks_bad) == 1, ranks_bad
    assert min(ranks_ok) == max(ranks_ok) == 1, ranks_ok

    def spread(build, npts):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            v = [endpoint(RadauIIA3Integrator, build, npts, off)
                 for off in (0.0, 0.31, 0.79)]
        return float(max(v) - min(v))

    ## (1) rank-constant: radau is at its own order, and refining collapses it
    s_lin = spread(constant_rank, 200)
    s_nl_200 = spread(nonlinear_constant_rank, 200)
    s_nl_400 = spread(nonlinear_constant_rank, 400)
    assert s_nl_200 < 1e-10, s_nl_200
    assert s_nl_400 < s_nl_200 / 10.0, (s_nl_200, s_nl_400)

    ## (2) the confound control: nonlinearity ALONE costs nothing
    assert s_nl_200 < 100.0 * max(s_lin, 1e-15), (s_lin, s_nl_200)

    ## (3) rank CHANGING: four orders worse, and the rate is gone
    s_bad_200 = spread(rank_changing, 200)
    s_bad_400 = spread(rank_changing, 400)
    assert s_bad_200 > 1e-6, s_bad_200
    assert s_bad_200 > 1e4 * s_nl_200, (s_nl_200, s_bad_200)
    assert s_bad_400 > s_bad_200 / 3.0, (s_bad_200, s_bad_400)


def test_one_netlist_returns_three_different_solutions_chosen_by_the_newton_seed():
    """⚠⚠ ONE NETLIST, ONE GRID, ONE TOLERANCE -- THREE DIFFERENT ANSWERS, and
    the choice is made SILENTLY by the Newton's initial guess.

    At a critical point of a DAE the solution can be genuinely non-unique
    (Lamour, März & Tischendorf Thm 3.53: "there are TWO solutions passing
    through").  `test_violating_the_im_D_hypothesis_...` next door measures a
    HARMLESS critical point -- order collapses, uniqueness survives.  This one
    is not harmless, and it needs TWO ingredients rather than one:

    * `rank C` DROPS -- `C = c0 V^2` vanishes at `V = 0`;
    * the equilibrium there is REPELLING -- a NEGATIVE conductance.

    With a passive conductance the equilibrium attracts, the field is one-sided
    Lipschitz, and forward uniqueness is safe however badly `C` degenerates.  A
    negative conductance is not exotic: it is what an oscillator's active
    device supplies.

    ⚠⚠ THE KNOB IS NOT THE GRID.  The analytic non-uniqueness becomes
    MULTIPLICITY OF ROOTS OF THE STEP EQUATION -- implicit Euler from
    `v_prev = 0` gives `z(c0 z^2/3h + g) = 0`, three roots when `g < 0` -- and
    the solver picks one silently, the SAME one on every grid.  A grid
    refinement or grid-offset probe returns a clean, convincing, wrong null.
    (Construction relayed from a peer session; measured here.)

    MEASURED, `V(t=1)` from `V = 0` by the first step's seed: -1.414744, 0, or
    +1.414744 at N = 3200, spread 2.83 that does NOT shrink (2.8405 / 2.8321 /
    2.8295 at N = 200/800/3200).  The non-trivial branches are exact --
    `w' = (3w)^(1/3)` integrates to `V(1) = sqrt(2) = 1.414214`.  With `g = +1`
    every seed returns identically 0.

    ⚠ A REPELLING EQUILIBRIUM CANNOT BE ARRIVED AT, which is why this starts
    ON it.  An earlier version drove an orbit "through" the point and found
    nothing -- on orbits where `min|V|` never fell below `|v0|`, because with
    `g < 0` the origin repels and a forward trajectory can only leave it.

    ⚠ AND THE STAGE PREDICTOR SHIPPED THE SAME DAY DOES NOT CHANGE THE CHOICE:
    it changed every Newton seed in this tree, and the seed is what selects the
    branch, so this is measured rather than assumed.  The first step from
    `V = 0` has no history, so the predictor declines and falls back.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from branch_selection import march

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ## (1) the repelling case: three distinct branches, and the spread
        ## does NOT shrink under refinement -- that is what says non-unique
        ## rather than inaccurate
        v200 = [march(-1.0, 200, s) for s in (-1.0, 0.0, 1.0)]
        v800 = [march(-1.0, 800, s) for s in (-1.0, 0.0, 1.0)]
        ## (2) the attracting control: uniqueness is safe, same degeneracy
        ctrl = [march(1.0, 200, s) for s in (-1.0, 0.0, 1.0)]
        ## (3) the predictor shipped today does not change the branch
        off = march(-1.0, 400, None, 'off', RadauIIA3Integrator)
        on = march(-1.0, 400, None, 'on', RadauIIA3Integrator)

    s200 = max(v200) - min(v200)
    s800 = max(v800) - min(v800)
    assert s200 > 2.0, v200
    assert s800 > 0.95 * s200, (s200, s800)          # does NOT shrink
    ## the non-trivial branches are +-sqrt(2), and the middle one is exactly 0
    assert abs(abs(v800[0]) - np.sqrt(2.0)) < 5e-3, v800
    assert abs(abs(v800[2]) - np.sqrt(2.0)) < 5e-3, v800
    assert v800[1] == 0.0, v800
    ## the control has the SAME vanishing C and is unique anyway
    assert max(ctrl) - min(ctrl) == 0.0, ctrl
    ## and today's predictor picks the same branch as the seed it replaced
    assert abs(on - off) < 1e-9, (off, on)


def test_the_branch_check_reports_a_multi_root_step_and_stays_quiet_otherwise():
    """The diagnostic for `test_one_netlist_returns_three_different_solutions_
    chosen_by_the_newton_seed`: `Transient.branch_check`, default ON.

    THE PROBLEM IT SOLVES.  Where `rank C` drops and the surviving dynamics
    repel, the step equation has several roots and the Newton takes whichever
    one its seed is nearest, silently.  There is no way to know a root was
    non-unique without looking for another one, so the check RE-SOLVES the same
    step from a perturbed seed -- which is why it needs a screen in front of it.

    ⚠⚠ THE OBVIOUS SCREEN WAS MEASURED AND REJECTED.  A peer session proposed
    firing on the STEP MATRIX ACQUIRING A NEGATIVE EIGENVALUE (`C/h + G < 0`
    needs `C` small AND `G` negative, so one number carries both ingredients)
    and verified it 6/6 against root counts on scalar and 2x2 systems.  On real
    MNA it fires on EVERYTHING -- measured, `min Re eig(J)` is negative on the
    index-2 C-V loop (-9.99e-07), the exponential fixture (-9.90e-07) and a van
    der Pol (-1.00e+00).  Two structural reasons: MNA WITH A VOLTAGE SOURCE IS
    A SADDLE-POINT SYSTEM, indefinite by construction, and AN OSCILLATOR'S `G`
    HAS A NEGATIVE EIGENVALUE BY DESIGN with no rank drop anywhere.  4 false
    fires out of 4 ordinary circuits.

    What is screened instead is the condition itself -- `rank C(x)` below the
    STRUCTURAL rank, which is "im D(t) is not time-invariant".  ⚠ It must be
    the structural rank and not a running maximum: on the degenerate branch
    `C` is identically zero for the whole run, so its rank never "drops", and
    screening against a running maximum reads QUIET on the very fixture that
    motivated this.

    ⚠ THE ASYMMETRY IS THE POINT.  The confirmation perturbs by a HEURISTIC
    magnitude, so it can MISS a second root; when it fires it has an actual
    second solution in hand.  A warning is evidence, silence is not.

    Asserted: it fires on the repelling fixture, is silent on the attracting
    control WHOSE `C` DEGENERATES IDENTICALLY (so the two-ingredient
    requirement is what separates them, not the rank drop alone), and does not
    fire on ordinary circuits -- ON EVERY SOLVE PATH, including the COUPLED
    one, which goes through `_stage_newton` rather than `_newton`.

    ⚠⚠ EXTENDING IT TO THE COUPLED PATH FOUND THREE DEFECTS, ALL MINE, AND THE
    LAST ONE IS THE INTERESTING ONE:

    * "fired" and "gave a direction" are different answers -- when `C`
      collapses ENTIRELY the screen returns `(True, None)`, and a
      `fired = direction` idiom reads that as "did not fire".  Zero screens on
      the fixture it was written for.
    * a FIXED-POINT test does not verify a root: if the solve hands its seed
      back, re-solving from that seed hands it back again and the test passes
      vacuously.  The BLOCK RESIDUAL is assembled and measured instead.
    * ⚠⚠ THE PERTURBATION WAS A GAUGE SHIFT.  The coupled path perturbs
      FULL-WIDTH stage vectors, and a direction of `ones` moves the REFERENCE
      NODE too -- a common-mode shift the circuit cannot see.  The solve leaves
      the pinned row alone, the "alternative" differs from the base only there,
      and its residual is EXACTLY ZERO because it is the same physical
      solution.  Measured: `alt` came back [0.7071, 0.7071] with `r_alt = 0.0`.
      FIVE FALSE ALARMS OUT OF FIVE on the attracting control, and a residual
      check could not catch it because the residual was genuinely zero.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.circuit import gnd as _gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (Gear2Integrator,
                                              RadauIIA3Integrator)
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from branch_selection import build
    from pycircuit.circuit.tests.test_stage_predictor import (_expg_fixture,
                                                             PER)

    def counts(g, npts=5, tend=1.0 / 40, cls=Gear2Integrator):
        cir = build(g)
        tr = Transient(cir, integrator=cls(), reltol=1e-12)
        tr.irefnode = cir.get_node_index(_gnd)
        x = np.zeros(cir.n)
        tr.epar.t = 0.0
        tr._begin_run(x, cir.n)
        h = tend / npts
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for j in range(1, npts + 1):
                tr._dt_last = tr._dt if j > 1 else None
                tr._dt = h
                tr.epar.t = j * h
                x, _f, _J, _ = tr.solve_timestep(x, j * h)
                tr._push_history(x)
        assert getattr(tr, '_branch_error', None) is None, tr._branch_error
        return (getattr(tr, 'branch_screens', 0),
                getattr(tr, 'branch_points', 0))

    ## (1)+(2), on EVERY solve path: the multistep one, the two
    ## DIRK-sequential ones, the GLM stage one, and the COUPLED one -- which
    ## goes through `_stage_newton` rather than `_newton` and so needed its own
    ## wiring.  The repelling fixture must fire; the attracting control, whose
    ## `C` degenerates IDENTICALLY, must not.  Running both on each path is
    ## what says the check tests MULTIPLICITY and not merely the rank drop.
    from pycircuit.circuit.integrator import (ESDIRK43Integrator,
                                              TRBDF2Integrator,
                                              GLM3Integrator)
    for cls, name in ((RadauIIA3Integrator, 'radau (coupled)'),
                      (ESDIRK43Integrator, 'esdirk43'),
                      (TRBDF2Integrator, 'trbdf2'),
                      (Gear2Integrator, 'gear'),
                      (GLM3Integrator, 'glm3')):
        s_bad, p_bad = counts(-1.0, cls=cls)
        assert s_bad > 0, (name, s_bad)
        assert p_bad > 0, (name, p_bad)
        s_ok, p_ok = counts(1.0, cls=cls)
        assert s_ok > 0, (name, s_ok)
        assert p_ok == 0, (name, 'FALSE ALARM', p_ok)

    ## (3) an ordinary circuit never even reaches the expensive half
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_expg_fixture(PER), integrator=Gear2Integrator(),
                       reltol=1e-10)
        tr.solve(refnode=_gnd, tend=PER, timestep=PER / 100,
                 fixed_timestep=True)
    ## (attributes, not `getattr(..., 0)`: that default read 0 while the
    ## check was switching itself off on every full solve, 2026-09-27)
    assert tr.statistics.branch_screens == 0
    assert tr.statistics.branch_points == 0
    assert getattr(tr, '_branch_error', None) is None, tr._branch_error

    ## (4) and it is switchable
    prev = Transient.branch_check
    try:
        Transient.branch_check = 'off'
        s_off, p_off = counts(-1.0)
    finally:
        Transient.branch_check = prev
    assert s_off == 0 and p_off == 0, (s_off, p_off)


def test_the_branch_check_reports_on_a_full_transient_solve(caplog):
    """`branch_check` on the path a user takes -- `Transient.solve` -- and
    not only on the hand-driven marches the test above counts on.

    ⚠⚠ UNTIL 2026-09-27 IT NEVER REPORTED THERE.  `TransientStatistics` has
    `__slots__`, and `branch_screens` / `branch_points` were not among them,
    so the first screen of every full solve raised AttributeError in
    `_branch_count`; the check's own except caught it, logged "the branch
    check itself failed" and switched the check off for the object's life.
    The marches have no `statistics` object and count on the instance, so
    every test passed, and the one full solve asserted
    `getattr(tr.statistics, 'branch_screens', 0) == 0` -- which the missing
    attribute satisfied.

    Asserted on the fixture of `test_one_netlist_returns_three_different_
    solutions_chosen_by_the_newton_seed`, for an LMM and the coupled Radau:
    the repelling equilibrium is REPORTED (screens, confirmed points, and the
    logged warning), the attracting control whose `C` degenerates identically
    is screened and NOT confirmed, and the check is still on afterwards.
    """
    import logging
    import os
    import sys
    from pycircuit.circuit.circuit import gnd as _gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (Gear2Integrator,
                                              RadauIIA3Integrator)
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from branch_selection import build

    for cls in (Gear2Integrator, RadauIIA3Integrator):
        for g in (-1.0, 1.0):
            caplog.clear()
            tr = Transient(build(g), integrator=cls(), reltol=1e-10)
            with caplog.at_level(logging.WARNING):
                tr.solve(refnode=_gnd, tend=1.0, timestep=1.0 / 50,
                         fixed_timestep=True)
            st = tr.statistics
            logged = [r.getMessage() for r in caplog.records]
            name = (cls.__name__, g)
            assert getattr(tr, '_branch_error', None) is None, (name, tr._branch_error)
            assert not any('branch check itself failed' in m for m in logged), name
            assert tr.branch_check == 'on', name
            assert st.branch_screens > 0, (name, st.branch_screens)
            if g < 0:
                assert st.branch_points > 0, (name, st.branch_points)
                assert any('MORE THAN ONE ROOT' in m for m in logged), (name, logged)
            else:
                assert st.branch_points == 0, (name, 'FALSE ALARM', st.branch_points)
                assert not any('MORE THAN ONE ROOT' in m for m in logged), name

    ## ⚠ A REAL JUNCTION IS NOT A RANK DROP.  With the check live, a
    ## half-wave rectifier screened 2086 of 2500 steps and re-solved each
    ## (3.8-4.9x the run on three junction circuits): the collapse scale came
    ## from random states in +-1 V, where the diode's diffusion capacitance
    ## is 0.217 F, so its real 1.7e-10 F junction read as zero.  The scale is
    ## the run's own now; a live junction is never screened.
    import pycircuit.circuit.elements_hdl as eh
    from pycircuit.circuit.elements import VSin
    rect = SubCircuit()
    a, b = rect.add_node('a'), rect.add_node('b')
    rect['V1'] = VSin(a, _gnd, va=5.0, freq=1e3)
    rect['D1'] = eh.DiodeSpiceHdl(a, b, IS=1.2e-14, rs=1.5, n=1.06, tt=4e-9,
                                  cjo=2.3e-12, vj=0.78, m=0.42, eg=1.11,
                                  xti=3.0, fc=0.5, bv=45.0, ibv=5e-6, kf=0.0,
                                  af=1.0, area=50.0, tnom=27.0)
    rect['Rl'] = R(b, _gnd, r=1e4)
    rect['Cl'] = C(b, _gnd, c=1e-6)
    rect.update_iparv()
    tr = Transient(rect)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(refnode=_gnd, tend=2e-3, timestep=2e-6)
    assert tr.statistics.accepted_steps > 100
    assert tr.statistics.branch_screens == 0, tr.statistics.branch_screens

    ## a failure of the check itself disables it for THAT run only: the next
    ## `solve` on the same object checks again (it used to overwrite
    ## `branch_check` with 'off' for the object's life)
    tr = Transient(build(-1.0), integrator=Gear2Integrator(), reltol=1e-10)

    def boom(*a, **kw):
        raise RuntimeError('screen failed')
    tr._branch_screen = boom
    with caplog.at_level(logging.WARNING):
        tr.solve(refnode=_gnd, tend=1.0, timestep=1.0 / 50, fixed_timestep=True)
    assert 'screen failed' in str(tr._branch_error)
    assert tr.statistics.branch_points == 0
    del tr._branch_screen
    tr.solve(refnode=_gnd, tend=1.0, timestep=1.0 / 50, fixed_timestep=True)
    assert tr._branch_error is None and tr.branch_check == 'on'
    assert tr.statistics.branch_points > 0


def test_the_branch_check_does_not_disturb_device_limiting_state():
    """⚠⚠ A DIAGNOSTIC THAT CHANGES THE SIMULATION IS A DEFECT, AND THIS ONE
    DID.  The confirmation re-solves a step from a perturbed seed, and every
    solve path in `transient.py` calls `cir.limit`, which for a junction device
    WRITES `_vlim` on the instance.  So the speculative solve left the limiting
    state at the ALTERNATIVE's value -- and `Diode.G` linearises around
    `_vlim`, so the NEXT step's Jacobian was then taken at the wrong point.

    MEASURED on a rank-dropping circuit carrying a diode: `_vlim` read 0.0 with
    `branch_check='off'` and 0.10166261963824502 with it on, while the step's
    own answer was unchanged.  A LATENT corruption -- it only bites once the
    screen fires, which is why the full suite never saw it.

    The fix puts the state back after every confirmation, whether or not it
    found anything -- from a snapshot, since 2026-09-24: the `limit(x, x)`
    re-sync it first used lands short above a junction's critical voltage
    (see the next test).  Asserted on both the single-Newton path and the
    coupled one, since they restore separately.
    """
    import warnings
    import numpy as np
    import sympy
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.circuit import SubCircuit, gnd as _gnd
    from pycircuit.circuit.elements import G as _G, Diode
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, ddt)
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (RadauIIA3Integrator,
                                              Gear2Integrator)

    _cc.default_toolkit = numeric

    class CubicCap(Behavioural):
        instparams = [Parameter(name='c0', desc='c', unit='F', default=1.0)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, ddt(c0 * b.V ** 3 / 3)),)  # noqa: F821

    def build():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['cq'] = CubicCap('a', _gnd, c0=1.0)
        c['g'] = _G('a', _gnd, g=-1.0)
        c['gd'] = _G('a', 'b', g=1e-3)
        c['d'] = Diode('b', _gnd, IS=1e-15)
        return c

    def march(cls, chk, npts=2, h=1.0 / 200):
        cir = build()
        tr = Transient(cir, integrator=cls(), reltol=1e-11)
        tr.branch_check = chk
        tr.irefnode = cir.get_node_index(_gnd)
        x = np.zeros(cir.n)
        tr.epar.t = 0.0
        tr._begin_run(x, cir.n)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for j in range(1, npts + 1):
                tr._dt_last = tr._dt if j > 1 else None
                tr._dt = h
                tr.epar.t = j * h
                x, _f, _J, _ = tr.solve_timestep(x, j * h)
                tr._push_history(x)
        return (np.asarray(x, dtype=float), getattr(cir['d'], '_vlim', None),
                getattr(tr, 'branch_points', 0))

    for cls, name in ((RadauIIA3Integrator, 'radau (coupled)'),
                      (Gear2Integrator, 'gear (single Newton)')):
        x_off, v_off, p_off = march(cls, 'off')
        x_on, v_on, p_on = march(cls, 'on')
        ## the check must have actually run, or this proves nothing
        assert p_on > 0 and p_off == 0, (name, p_on, p_off)
        ## and it must have left no trace
        assert v_off is not None and v_on is not None, (name, v_off, v_on)
        assert abs(float(v_on) - float(v_off)) < 1e-12, (name, v_off, v_on)
        assert np.max(np.abs(x_on - x_off)) < 1e-12, (name, x_on, x_off)


def test_the_branch_check_confirms_on_every_solve_path():
    """⚠ THE SILENT HOLE IS THE PROBLEM, NOT THE MISSING CONFIRMATION -- and
    the confirmation is no longer missing (2026-09-24).  Three paths --
    Radau's opt-in transform, the coupled PCNR step and the multistep PCNR
    step -- carried their own Newton, not re-enterable from a seed, so they
    only SCREENED (`branch_screens_unconfirmed`), and could not tell the
    repelling fixture from the attracting one.  And a fourth ran NOTHING: a
    DIRK or GLM stage solved by PCNR (trbdf2, esdirk43, the GLMs under
    `pcnr=True`) -- neither the screen nor the confirmation.

    The step equation is the same whichever Newton solves it, so each path
    now confirms with the limiting one: the multistep and stage paths
    through `_branch_after_solve` on their own residual, the coupled paths
    through the dense coupled Newton (`_coupled_stage_solver`), built only
    once the screen fires.  Measured on the repelling (g = -1) and
    attracting (g = +1) fixtures: every path confirms the first and screens
    but does not confirm the second, and none leaves an unconfirmed screen.
    """
    import os
    import sys
    import warnings
    import numpy as np
    from pycircuit.circuit.circuit import gnd as _gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit import integrator as _I
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from branch_selection import build
    from pycircuit.circuit.tests.test_stage_predictor import (_expg_fixture,
                                                              PER)

    def march(g, cls, transform=False, pcnr=False):
        cir = build(g)
        tr = Transient(cir, integrator=cls(), reltol=1e-12, pcnr=pcnr)
        tr._radau_use_transform = transform
        tr.irefnode = cir.get_node_index(_gnd)
        x = np.zeros(cir.n)
        tr.epar.t = 0.0
        tr._begin_run(x, cir.n)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for j in range(1, 4):
                tr._dt_last = tr._dt if j > 1 else None
                tr._dt = 1.0 / 200
                tr.epar.t = j / 200.0
                x, _f, _J, _ = tr.solve_timestep(x, j / 200.0)
                tr._push_history(x)
        return (getattr(tr, 'branch_screens', 0), getattr(tr, 'branch_points', 0),
                getattr(tr, 'branch_screens_unconfirmed', 0),
                getattr(tr, '_branch_error', None))

    paths = (('transform', _I.RadauIIA3Integrator, True, False),
             ('coupled PCNR', _I.RadauIIA3Integrator, False, True),
             ('multistep PCNR', _I.Gear2Integrator, False, True),
             ('DIRK PCNR', _I.TRBDF2Integrator, False, True),
             ('GLM PCNR', _I.GLM2Integrator, False, True))
    for name, cls, tf, pcnr in paths:
        s_bad, p_bad, u_bad, e_bad = march(-1.0, cls, tf, pcnr)
        assert e_bad is None, (name, e_bad)
        assert s_bad > 0 and p_bad > 0 and u_bad == 0, (name, s_bad, p_bad, u_bad)
        s_ok, p_ok, u_ok, _e = march(1.0, cls, tf, pcnr)
        assert s_ok > 0 and p_ok == 0 and u_ok == 0, (name, s_ok, p_ok, u_ok)

    ## and an ordinary circuit is quiet on both transform settings
    for tf in (True, False):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr = Transient(_expg_fixture(PER), integrator=_I.RadauIIA3Integrator(),
                           reltol=1e-10)
            tr._radau_use_transform = tf
            tr.solve(refnode=_gnd, tend=PER, timestep=PER / 100,
                     fixed_timestep=True)
        assert getattr(tr.statistics, 'branch_points', 0) == 0


def test_the_one_over_h_defect_amplification_is_LOCAL_not_propagated():
    """⚠⚠ THE QUESTION THAT WAS OPEN ALL SESSION, ANSWERED: the amplification
    is LOCAL.

    The filtered local-error estimate reads one order below its declared
    `EMBEDDED_ORDER + 1` on a DAE, because `J = C + a·h·G` is singular in `C`
    so `J⁻¹` behaves like `1/h` on the algebraic subspace. Whether that matters
    depends entirely on whether the amplified part PROPAGATES. Lamour, März &
    Tischendorf §8.4 note (6) says it does not, at index ≤ 2 with a properly
    stated leading term. This measures it here rather than citing it.

    ⚠⚠ TWO THINGS HAD TO BE RIGHT AT ONCE, and fixing either alone still reads
    nothing — which is why three earlier attempts failed:

    * the defect must go in a direction `‖J⁻¹‖` ACTUALLY AMPLIFIES, the left
      singular vector of `J` for `σ_min`. A defect in that operator's
      nullspace is amplified by exactly nothing, and `|e| = δ` is then the
      CORRECT answer — which is what the earlier attempts were measuring;
    * and it must be probed BELOW the `‖C‖/‖G‖` turn. Above it the reactive
      term is negligible, `J` is effectively resistive, and there is no
      amplification anywhere to find.

    MEASURED, `δ = 1e-9` injected at ONE step through `provided_function` (a
    defect in the residual — the theorem's `q_ni`), differenced against the
    undisturbed run of the SAME discretisation so truncation cancels exactly:

    ===========  =====================  ==========================
    fixture      amplification at n0    tail after μ steps
    ===========  =====================  ==========================
    C-V loop     −1.00, −1.00 (μ=2)     +2.00, +1.98 → shrinks h²
    ExpG         −0.00, −0.00 (μ=1)     +1.00, +1.00 → shrinks h
    ===========  =====================  ==========================

    The INSTRUMENT-ALIVE half is the amplification column: the index-2
    fixture amplifies as `δ/h` to two decimals, exactly Prop 8.10's
    `h^-(μ-1)`, and the index-1 fixture is correctly flat. Without that, a
    small tail proves nothing.

    ⚠ The arm of note (6) covering a variable-coefficient nonlinear MNA is the
    index-2 one, and the margin to the PROPAGATING index-3 case is one index
    level. Nothing here tests index 3.
    """
    import os
    import sys
    import warnings
    import numpy as np
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', '..', '..', 'benchmarks'))
    from defect_locality import locality_below_the_turn, index1, index2

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        amp2, tail2 = locality_below_the_turn(index2, 'index-2', 2)
        amp1, tail1 = locality_below_the_turn(index1, 'index-1', 1)

    ## (1) INSTRUMENT ALIVE: index 2 amplifies as delta/h, index 1 does not
    for v in amp2:
        assert abs(v - (-1.0)) < 0.1, ('index-2 amplification', amp2)
    for v in amp1:
        assert abs(v) < 0.1, ('index-1 amplification', amp1)

    ## (2) LOCAL: what survives mu steps SHRINKS with refinement.  A
    ## propagated term would keep the amplification's own scaling (negative).
    for v in tail2:
        assert v > 0.5, ('index-2 tail does not decay -- PROPAGATED?', tail2)
    for v in tail1:
        assert v > 0.5, ('index-1 tail does not decay', tail1)


def test_the_branch_check_solves_the_step_equation_and_restores_device_state_exactly():
    """Two defects of the branch check's confirmation, on a junction ABOVE
    its critical voltage -- which the test above cannot reach: its diode sits
    at 0 V, where `pnjlim` never clamps.

    (1) THE SPECULATIVE SOLVE PASSED NO LIMITER, so no `cir.limit` ran and a
    stateful device (`Diode`, whose `i` and `G` read `_vlim`) stayed
    linearised at the ACCEPTED point: the solve and its root test both saw
    that linearisation.  Measured with a diode ON the rank-dropping node at
    0.85 V: the reported alternative, a = 1.0536, had residual 0 in the
    linearised system and 489 A in the step equation.  With the step's own
    limiter it is a = 0.86584, a root of the step equation.

    (2) THE RESTORE WAS `limit(x, x)`, which clamps against the STORED
    state: from a speculative 1.054 V it landed at 0.790, not 0.85.  It is a
    snapshot now.  Asserted by leaving the state far away inside the
    confirmation, as a speculative Newton that wandered would, and requiring
    the step to end with the state it has with the check off.
    """
    import types
    import warnings
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.circuit import SubCircuit, gnd as _gnd
    from pycircuit.circuit.elements import G as _G, Diode, VS as _VS, IS as _IS
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, ddt)
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (RadauIIA3Integrator,
                                              Gear2Integrator,
                                              TrapezoidalIntegrator)

    _cc.default_toolkit = numeric
    v0, IS_, k0 = 0.85, 1e-15, 10.0
    VT = 1.380649e-23 * 300.15 / 1.602176634e-19

    class CubicCapAt(Behavioural):
        instparams = [Parameter(name='c0', desc='c', unit='F', default=1.0),
                      Parameter(name='v0', desc='v0', unit='V', default=0.0)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, ddt(c0 * (b.V - v0) ** 3 / 3)),)  # noqa: F821

    def build():
        ## `a` has C = 0 at v0 and a net NEGATIVE conductance there (k0 minus
        ## the diode's own), with the diode's current at v0 fed in: an
        ## equilibrium on a rank drop, so the step equation has several roots
        c = SubCircuit()
        c.add_node('a')
        c.add_node('r')
        c['cq'] = CubicCapAt('a', _gnd, c0=1.0, v0=v0)
        c['d'] = Diode('a', _gnd, IS=IS_)
        c['vr'] = _VS('r', _gnd, v=v0)
        c['gk'] = _G('a', 'r', g=-k0)
        c['i0'] = _IS(_gnd, 'a', i=IS_ * (np.exp(v0 / VT) - 1.0))
        return c

    def march(cls, chk, wander=None, npts=3, h=1.0 / 200):
        cir = build()
        d = cir['d']
        tr = Transient(cir, integrator=cls(), reltol=1e-11)
        tr.branch_check = chk
        tr.irefnode = cir.get_node_index(_gnd)
        alts = []
        confirm, scan = tr._branch_confirm, tr._branch_coupled_scan

        def spy_confirm(self, func, x_res, direction):
            alt = confirm(func, x_res, direction)
            if alt is not None:
                ## the step equation ITSELF at `alt`: the diode evaluated at
                ## the point (no stored state), then put back
                saved = dict(d.__dict__)
                d.__dict__.pop('_vlim', None)
                F, _J = func(alt)
                d.__dict__.clear()
                d.__dict__.update(saved)
                alts.append((float(np.asarray(alt, dtype=float)[0]),
                             float(np.max(np.abs(np.asarray(F, dtype=float))))))
            if wander is not None:
                d.__dict__['_vlim'] = wander
            return alt

        def spy_scan(self, *a, **k):
            gap = scan(*a, **k)
            if wander is not None:
                d.__dict__['_vlim'] = wander
            return gap
        tr._branch_confirm = types.MethodType(spy_confirm, tr)
        tr._branch_coupled_scan = types.MethodType(spy_scan, tr)
        x = np.zeros(cir.n)
        x[cir.get_node_index('a')] = v0
        x[cir.get_node_index('r')] = v0
        tr.epar.t = 0.0
        tr._begin_run(x, cir.n)
        vl = []
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for j in range(1, npts + 1):
                tr._dt_last = tr._dt if j > 1 else None
                tr._dt = h
                tr.epar.t = j * h
                x, _f, _J, _ = tr.solve_timestep(x, j * h)
                tr._push_history(x)
                vl.append(float(d.__dict__['_vlim']))
        return (np.asarray(x, dtype=float), vl, alts,
                getattr(tr, 'branch_points', 0), getattr(tr, 'branch_screens', 0))

    for cls in (Gear2Integrator, TrapezoidalIntegrator):
        name = cls.__name__
        x_off, v_off, _a, p_off, _s = march(cls, 'off')
        x_on, v_on, alts, p_on, _s = march(cls, 'on')
        ## the check ran, and fired, or this proves nothing
        assert p_on > 0 and p_off == 0 and alts, (name, p_on, p_off, alts)
        ## (1) every alternative it reports is a root of the step equation
        for a, F in alts:
            assert abs(a - v0) > 1e-6 and F < 1e-9, (name, a, F)
        assert np.max(np.abs(x_on - x_off)) < 1e-12, (name, x_on, x_off)
        assert v_on == v_off, (name, v_on, v_off)

    ## (2) whatever the speculative solve leaves, the state goes back, on
    ## both solve paths (the coupled one is radau's).  1.0536 is where the
    ## unlimited solve's "alternative" sat, 0.2 V beyond the clamp's reach of
    ## 0.85; 0 V is a junction a wandering Newton turned off.  (The coupled
    ## path's old restore re-synced once PER STAGE, and three calls walk
    ## back from 1.0536 -- 0.790, 0.821, 0.850 -- but not from 0 V.)
    for cls in (Gear2Integrator, RadauIIA3Integrator):
        name = cls.__name__
        x_off, v_off, _a, _p, _s = march(cls, 'off')
        for wander in (1.0536, 0.0):
            x_on, v_on, _a, _p, screens = march(cls, 'on', wander=wander)
            assert screens > 0, (name, screens)
            assert v_on == v_off, (name, wander, v_on, v_off)
            assert np.max(np.abs(x_on - x_off)) < 1e-12, (name, wander, x_on, x_off)


def test_the_transient_lands_declared_state_events_and_its_period_stops_jittering():
    """E7 (2026-09-22, Andreas: "Add fixing the transient at events to our
    list"): `Transient` cuts an accepted step to a declared crossing
    (`Circuit.state_events()`, a VSwitch's two window edges) by a secant
    on the fraction and restarts the history there, as at a source corner.

    ⚠ WHAT IT BUYS, MEASURED on the sharp comparator oscillator against the
    windowed exact period (25 periods): with VSwitch's compact transition
    the LTE controller already localises the crossing to the tolerance,
    so an LMM's MEAN period error is its own and landing does not change
    it -- gear +5.4e-3 / 1.4e-3 / 3.7e-4 / 9.1e-5 unlanded, +6.7e-3 /
    1.8e-3 / 4.2e-4 / 9.6e-5 landed at reltol 1e-4 .. 1e-7, restart or
    not.  What landing removes is the period-to-period JITTER, the "where
    in the step" randomness: gear's spread 1.3e-4 -> 1.5e-5 at 1e-6
    (6.5e-3 -> 3e-4 at 1e-4); trap's -1.1e-5 mean was a cancellation
    inside a 1.2e-3 spread and becomes +1.9e-4 +- 1.3e-5.  Cost: ~2 secant
    cuts per landing, +30 % steps.  Radau, one-step, is at -1.7e-8
    unlanded; wired into its `_run_rk_adaptive` loop the same day, landing
    trims its spread 4.8e-6 -> 8.2e-7 (1e-4) and 3.7e-8 -> 1.9e-9 (1e-6),
    mean unchanged, steps +40-130 %; trbdf2 gains nothing in spread.
    Pinned: four landings per period, each on a window edge to 5 % of
    the window (1e-3 of the step; measured 1.1 %), the landed spread
    below 5e-5 and the unlanded above it, `state_events=False` landing
    nothing; no landing in the first step (an initial condition is not a
    solved state); the Runge-Kutta loop, wired the same day, landing the
    same edges for radau; and the coupled-LTE loop, wired last, where
    landing does buy accuracy (gear coupled: mean +1.4e-3 -> +1.9e-4,
    spread 2.7e-3 -> 2.7e-6 at reltol 1e-6)."""
    import warnings as _w
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    T_ex = 1.3918372887e-6
    cir = _comparator_relaxation_oscillator()
    names = [str(n_) for n_ in cir.nodes]
    x0 = np.zeros(cir.n)
    x0[names.index('c')] = 1.0
    ifb1, iref = names.index('fb1'), names.index('ref')
    out = {}
    for se in (True, False):
        tr = Transient(cir, toolkit=circuit.numeric, reltol=1e-6, state_events=se)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=14 * T_ex, timestep=T_ex / 200, x0=x0)
        tt = np.asarray(res.sweep_values, dtype=float)
        xx = np.asarray(res.x, dtype=float)
        d = xx[ifb1] - xx[iref]
        up = np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))
        tc = tt[up] - d[up] * (tt[up + 1] - tt[up]) / (d[up + 1] - d[up])
        per = np.diff(tc)[-6:]
        out[se] = (np.ptp(per) / T_ex, np.mean(per) / T_ex - 1.0, tr, tt, d)
    spread_l, mean_l, tr_l, tt_l, d_l = out[True]
    spread_u, mean_u, tr_u, _tt, _d = out[False]
    assert spread_l < 5e-5 and spread_u > 5e-5, (spread_l, spread_u)
    assert abs(mean_l) < 2e-3 and abs(mean_u) < 2e-3, (mean_l, mean_u)
    assert tr_u.statistics.state_events_hit == 0 and not tr_u.event_times
    ev = np.asarray(tr_l.event_times, dtype=float)
    ## four edges per period once the orbit has settled (the first period starts off it)
    settled = ev[ev > 4 * T_ex]
    assert 4 * 9 <= len(settled) <= 4 * 10 + 2, len(settled)
    assert tr_l.statistics.state_events_hit == len(ev)
    ## every landed time is a grid point ON a window edge: |fb1 - ref| = 1e-4
    ## to 5 % of the window (EVENT_LAND_RTOL is 1e-3 of the STEP, a fifth of
    ## this 0.04 ns window at a 7 ns step; measured 1.1 %)
    for te in settled:
        j = int(np.argmin(np.abs(tt_l - te)))
        assert abs(tt_l[j] - te) < 1e-12 * T_ex
        assert abs(abs(d_l[j]) - 1e-4) < 2e-4 * 5e-2, (te, d_l[j])
    ## no landing in the first step: the initial condition is not a solved state
    assert ev[0] > 0.5 * T_ex, ev[:3]
    ## the Runge-Kutta loop, wired the same day: radau at reltol 1e-4 lands
    ## the same four edges (to 1e-5 of the window), spread 2.0e-6 -> 3.9e-7
    ## over these 14 periods (4.8e-6 -> 8.2e-7 over 25), mean unchanged
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    rk = {}
    for se in (True, False):
        tr = Transient(cir, toolkit=circuit.numeric, reltol=1e-4, state_events=se,
                       integrator=RadauIIA3Integrator())
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=14 * T_ex, timestep=T_ex / 200, x0=x0)
        tt = np.asarray(res.sweep_values, dtype=float)
        xx = np.asarray(res.x, dtype=float)
        d = xx[ifb1] - xx[iref]
        up = np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))
        tc = tt[up] - d[up] * (tt[up + 1] - tt[up]) / (d[up + 1] - d[up])
        per = np.diff(tc)[-6:]
        rk[se] = (np.ptp(per) / T_ex, np.mean(per) / T_ex - 1.0, tr, tt, d)
    ## measured over the last 6 of 14 periods: 3.9e-7 landed, 2.0e-6 unlanded
    assert rk[True][0] < 1e-6 and rk[False][0] > 1e-6, (rk[True][0], rk[False][0])
    assert abs(rk[True][1]) < 1e-5 and abs(rk[False][1]) < 1e-5
    assert rk[False][2].statistics.state_events_hit == 0
    ev_rk = np.asarray(rk[True][2].event_times, dtype=float)
    settled_rk = ev_rk[ev_rk > 4 * T_ex]
    assert 4 * 9 <= len(settled_rk) <= 4 * 10 + 2 and ev_rk[0] > 0.5 * T_ex, (len(settled_rk), ev_rk[:2])
    tt_rk, d_rk = rk[True][3], rk[True][4]
    for te in settled_rk:
        j = int(np.argmin(np.abs(tt_rk - te)))
        assert abs(tt_rk[j] - te) < 1e-12 * T_ex
        assert abs(abs(d_rk[j]) - 1e-4) < 2e-4 * 5e-2, (te, d_rk[j])
    ## the coupled-LTE loop (`coupled_lte=True`), wired last -- and the one
    ## place landing buys ACCURACY, not only consistency: the coupled (x, h)
    ## solve's own step unknown walked through the crossing badly.  gear,
    ## reltol 1e-6, 14 periods: mean +1.4e-3 -> +1.9e-4, spread 2.7e-3 ->
    ## 2.7e-6 (at 1e-4: +2.1e-2 -> +3.6e-3, 2.4e-2 -> 3.3e-4); the handed
    ## step is HELD during a cut (`hold_h`), edges to 1-4 % of the window.
    cp = {}
    for se in (True, False):
        tr = Transient(cir, toolkit=circuit.numeric, reltol=1e-6, state_events=se)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=14 * T_ex, timestep=T_ex / 200, x0=x0, coupled_lte=True)
        tt = np.asarray(res.sweep_values, dtype=float)
        xx = np.asarray(res.x, dtype=float)
        d = xx[ifb1] - xx[iref]
        up = np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))
        tc = tt[up] - d[up] * (tt[up + 1] - tt[up]) / (d[up + 1] - d[up])
        per = np.diff(tc)[-6:]
        cp[se] = (np.ptp(per) / T_ex, np.mean(per) / T_ex - 1.0, tr, tt, d)
    assert cp[True][0] < 1e-5 and cp[False][0] > 1e-4, (cp[True][0], cp[False][0])
    assert abs(cp[True][1]) < 5e-4 and abs(cp[False][1]) > 5e-4, (cp[True][1], cp[False][1])
    assert cp[False][2].statistics.state_events_hit == 0
    ev_cp = np.asarray(cp[True][2].event_times, dtype=float)
    settled_cp = ev_cp[ev_cp > 4 * T_ex]
    assert 4 * 9 <= len(settled_cp) <= 4 * 10 + 2 and ev_cp[0] > 0.5 * T_ex, (len(settled_cp), ev_cp[:2])
    tt_cp, d_cp = cp[True][3], cp[True][4]
    for te in settled_cp:
        j = int(np.argmin(np.abs(tt_cp - te)))
        assert abs(tt_cp[j] - te) < 1e-12 * T_ex
        assert abs(abs(d_cp[j]) - 1e-4) < 2e-4 * 5e-2, (te, d_cp[j])


def test_the_runge_kutta_transient_loop_lands_source_corners():
    """Item 2 of the 2026-09-22 list: `_run_rk_adaptive` (radau, trbdf2,
    esdirk in the transient) stepped OVER a source's corners -- it had no
    breakpoint handling at all, where `_solve` cuts its step to the next
    `next_event`.  Now it lands them the same way.  Measured against the
    RC's EXACT response to the VPulse's ramps (a 0.1 ns edge into tau =
    1 us): radau unlanded 7.1e-5 / 5.1e-7 at reltol 1e-4 / 1e-6, landed
    1.3e-8 / 8.1e-10 with FEWER steps (109 -> 79, 144 -> 106); trbdf2
    gains nothing (2.0e-4 -> 2.1e-4, 1.1e-5 -> 1.1e-5: its own second
    order); gear, whose loop always landed, would read 8.9e-3 / 5.9e-4
    without.  Pinned: radau at reltol 1e-4 below 1e-7 with all four
    corners hit, and above 1e-5 with `next_event` silenced."""
    import warnings as _w
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    circuit.default_toolkit = circuit.numeric
    TD, TR, PW, TF, TAU = 2e-6, 1e-10, 5e-6, 1e-10, 1e-6

    def build():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=TD, tr=TR, tf=TF, pw=PW, per=20e-6)
        c['r'] = R('a', 'b', r=1e3)
        c['c'] = C('b', gnd, c=1e-9)
        return c

    def ramp(t, t0, w):
        s = np.clip(t - t0, 0, None)
        return (s - TAU * (1 - np.exp(-s / TAU))) / w

    def exact(t):
        t = np.asarray(t, dtype=float)
        return (ramp(t, TD, TR) - ramp(t, TD + TR, TR)
                - ramp(t, TD + TR + PW, TF) + ramp(t, TD + TR + PW + TF, TF))

    errs = {}
    for landed in (True, False):
        c = build()
        names = [str(n_) for n_ in c.nodes]
        if not landed:
            c.next_event = lambda t: float('inf')
        tr = Transient(c, toolkit=circuit.numeric, reltol=1e-4, integrator=RadauIIA3Integrator())
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=10e-6, timestep=0.2e-6)
        tt = np.asarray(res.sweep_values, dtype=float)
        vb = np.asarray(res.x, dtype=float)[names.index('b')]
        errs[landed] = float(np.max(np.abs(vb - exact(tt))))
        if landed:
            assert tr.statistics.breakpoints_hit == 4, tr.statistics.breakpoints_hit
            for tc in (TD, TD + TR, TD + TR + PW, TD + TR + PW + TF):
                assert np.min(np.abs(tt - tc)) < 1e-15, tc
        else:
            assert tr.statistics.breakpoints_hit == 0
    assert errs[True] < 1e-7 and errs[False] > 1e-5, errs
