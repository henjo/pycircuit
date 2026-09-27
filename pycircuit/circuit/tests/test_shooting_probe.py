"""Shooting tests: shooting probe.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_elements import _PllMultPd
from pycircuit.circuit.tests._shooting_fixtures import (_injection_lock_edge,
    _vdp_injected)


@pytest.mark.slow
def test_probe_shooting_finds_the_orbit_and_screens_for_instability():
    """B5 built: Bizzarri's probe, and the 2x2 power-flow screen.

    A periodic voltage source across a node, with `(A, f)` solved so the
    probe's OWN fundamental current vanishes -- at which point it sources
    nothing and can be removed.  The probe makes the circuit NON-AUTONOMOUS,
    so the period is known, there is no phase condition, and no `T = 0`
    trivial root to fall into.

    ⚠ IT IS NOT A CONVERGENCE AID.  The paper's own flagship high-Q Pierce
    example was solved with "conventional SH" and a tentative inductor
    current, not the probe.  What it buys is the sweep -- unstable cycles,
    coexisting solutions, a stability screen.

    ⚠⚠ WHAT IT COSTS, MEASURED, AND IT IS INHERENT RATHER THAN A DEFECT.
    A single-tone probe forces a SINUSOID.  A non-sinusoidal orbit cannot
    null the probe's whole current, only its FUNDAMENTAL, so what this
    returns is the first-harmonic-balance (describing-function) solution.
    On van der Pol the error is clean and QUADRATIC in harmonic content::

        mu     autonomous f   probe f     df/f        THD      df/f / THD^2
        0.10   0.159053       0.159134    +5.06e-04   0.0112       4.05
        0.30   0.158304       0.159134    +5.24e-03   0.0361       4.02
        1.00   0.150229       0.159134    +5.93e-02   0.1192       4.17

    ⚠ AND THE PROBE FREQUENCY IS THE SAME AT EVERY `mu` -- 0.159134, which
    is the LC resonance `1/(2 pi sqrt(LC))`.  That is not a bug either:
    `mu (u - u^3/3)` is odd and memoryless, so its describing function is
    purely REAL and shifts no phase, and first-harmonic balance therefore
    MUST land on the linear resonance.  The true frequency moves away from
    it as harmonics grow.  A gate that only checked "the probe converged"
    would have accepted a 5.9% frequency error at `mu = 1` without noticing.

    ⚠⚠ PROBE PLACEMENT IS CIRCUIT-SPECIFIC, AND ITS FAILURE IS NOT A SOLVER
    FAILURE.  Across van der Pol's only node with no series resistance, the
    inductor's DC current is unconstrained once `v` is forced, so a whole
    family satisfies periodicity and the shooting Jacobian is SINGULAR.
    Measured: periodicity error 2.11e-15 -- already a periodic solution --
    reported as `converged = False`.  `degenerate_placement` names that
    pairing instead of leaving a correct answer labelled non-convergent.
    """
    import warnings as _w
    from pycircuit.circuit.shooting import ProbeShooting
    circuit.default_toolkit = circuit.numeric

    def vdp(mu, rs):
        def build():
            c = SubCircuit()
            c.add_node('v')
            c['C'] = C('v', gnd, c=1.0)
            if rs > 0:
                c.add_node('x')
                c['RL'] = R('v', 'x', r=rs)
                c['L'] = L('x', gnd, L=1.0)
            else:
                c['L'] = L('v', gnd, L=1.0)
            c['B'] = BSource('v', gnd, gnd, 'v',
                             i_func=lambda u: mu * (u - u ** 3 / 3.0))
            return c
        return build

    f_lc = 1.0 / (2.0 * np.pi)

    ## (1) THE DEGENERATE PLACEMENT, detected as such.
    bad = ProbeShooting(vdp(1.0, 0.0), 'v', npts=200)
    deg, perr, conv = bad.degenerate_placement(2.0, f_lc)
    assert deg and not conv and perr < 1e-10, \
        'a placement leaving the inductor DC free should read degenerate ' \
        '(got degenerate=%r periodicity=%.2e converged=%r)' % (deg, perr, conv)
    ok = ProbeShooting(vdp(1.0, 1e-2), 'v', npts=200)
    deg2, _p2, conv2 = ok.degenerate_placement(2.0, f_lc)
    assert conv2 and not deg2, \
        'a series resistance removes the free mode, so this must converge'

    ## (2) NEARLY SINUSOIDAL: the probe must agree with the AUTONOMOUS solve,
    ## which is the only correctness evidence here -- converging to something
    ## proves nothing.
    mu = 0.1
    ref = PSS(vdp(mu, 1e-2)(), method='gear', reltol=1e-11)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ref.solve(period=6.6634, timestep=6.6634 / 300,
                  x0=np.array([2.0, 0.0, 0.0]), maxiterations=60)
    assert ref.converged, 'the autonomous reference did not converge'
    Xa = np.asarray(ref.waveform[1], dtype=float)
    f_ref = 1.0 / ref.period
    amp_ref = 0.5 * (Xa[0].max() - Xa[0].min())

    ps = ProbeShooting(vdp(mu, 1e-2), 'v', npts=300)
    A, f, info = ps.solve(amp_ref, f_ref, tol=1e-8, maxiter=12)
    assert info['converged'], 'the probe solve did not converge: %r' % (info,)
    assert abs(A - amp_ref) / amp_ref < 5e-3, \
        'probe amplitude %.6f against autonomous %.6f' % (A, amp_ref)
    assert abs(f - f_ref) / f_ref < 5e-3, \
        'probe frequency %.6f against autonomous %.6f' % (f, f_ref)

    ## (3) ⚠ AND THE FIRST-HARMONIC LIMIT IS ASSERTED, not left implicit: the
    ## probe lands on the LC resonance because the nonlinearity is odd and
    ## memoryless.  If this ever stops holding, the accuracy law above is
    ## wrong and the docstring must be re-measured.
    assert abs(f - f_lc) / f_lc < 2e-3, \
        'the probe should sit on the LC resonance %.6f for an odd memoryless ' \
        'nonlinearity, got %.6f' % (f_lc, f)

    ## (4) THE POWER-FLOW SCREEN on a circuit known stable.  ⚠ ONE-DIRECTIONAL:
    ## only `P > 0 => unstable` is proven, so a non-positive P means NOT
    ## DETECTED, never "stable".  Asserting the reverse would be asserting
    ## something the authors explicitly say is unproven.
    P, pinfo = ps.power_flow(A, f)
    assert not pinfo['unstable'], \
        'the power-flow screen flagged a stable van der Pol as unstable ' \
        '(P = %.6e)' % P
    assert pinfo['symmetric_part'].shape == (2, 2)
    assert np.allclose(pinfo['symmetric_part'], pinfo['symmetric_part'].T), \
        'the screen contracts the SYMMETRIC part, so it must be symmetric'


@pytest.mark.slow
def test_even_harmonic_pruning_must_be_measured_and_never_assumed():
    """⚠⚠ THE CHEAP OPTIMISATION THAT SILENTLY RETURNS THE WRONG ORBIT.

    A multi-tone probe costs `2K+1` PSS solves an iteration, and on van der
    Pol the even tones look like pure waste: `K=2` returns `A_2 = 2.3e-13` and
    a frequency identical to `K=1` in every printed digit.  It is tempting to
    drop them.

    **THAT HOLDS ONLY FOR A HALF-WAVE SYMMETRIC CIRCUIT.**  Add an even term
    `beta u^2` to the same nonlinearity and the even content is real::

        beta    H2/H1       H3/H1       H4/H1
        0.00    7.080e-16   1.168e-01   2.710e-16   <- symmetric
        0.05    2.649e-02   1.162e-01   9.239e-03
        0.20    1.058e-01   1.068e-01   3.600e-02   <- H2 EQUALS H3
        0.50    2.613e-01   5.994e-02   7.610e-02   <- H2 is 4x H3

    ⚠⚠ AND THE FAILURE IS SILENT.  On the asymmetric circuit (autonomous
    `f = 0.148220`) BOTH tone sets converge::

        tones       f           df/f        converged
        [1, 2, 3]   0.148753    +3.60e-03   True
        [1, 3]      0.150172    +1.32e-02   True      <- 3.7x worse

    ⚠⚠⚠ AND `0.150172` IS THE SYMMETRIC CIRCUIT'S OWN ANSWER, to every
    printed digit.  Dropping the even tones does not merely lose accuracy --
    it makes the probe STRUCTURALLY BLIND to `beta`, so it returns the orbit
    of a different circuit and reports convergence.  That is why
    `even_harmonic_content` measures instead of assuming, and why `tones` has
    no clever default.
    """
    import warnings as _w
    from pycircuit.circuit.shooting import ProbeShooting
    circuit.default_toolkit = circuit.numeric
    RS = 1e-2

    def fac(beta):
        def build():
            c = SubCircuit()
            c.add_node('v')
            c.add_node('x')
            c['C'] = C('v', gnd, c=1.0)
            c['RL'] = R('v', 'x', r=RS)
            c['L'] = L('x', gnd, L=1.0)
            c['B'] = BSource('v', gnd, gnd, 'v',
                             i_func=lambda u, b=beta: (u - u ** 3 / 3.0) + b * u ** 2)
            return c
        return build

    f_lc = 1.0 / (2.0 * np.pi)

    ## (1) THE MEASUREMENT MUST SEPARATE THE CASES, or the design is unusable.
    sym, _m = ProbeShooting(fac(0.0), 'v', npts=200).even_harmonic_content(f_lc)
    asym, _m2 = ProbeShooting(fac(0.2), 'v', npts=200).even_harmonic_content(f_lc)
    assert sym < 1e-10, \
        'the symmetric circuit should show no even content, got %.3e' % sym
    assert asym > 1e-2, \
        'the asymmetric circuit must show even content or the falsifier below ' \
        'proves nothing, got %.3e' % asym
    assert asym / max(sym, 1e-300) > 1e6, \
        'the measure must SEPARATE the cases, not merely order them'

    ## (2) THE SILENT FAILURE, asserted: pruning converges to a worse answer.
    ref = PSS(fac(0.2)(), method='gear', reltol=1e-11)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ref.solve(period=6.6634, timestep=6.6634 / 300,
                  x0=np.array([2.0, 0.0, 0.0]), maxiterations=60)
    assert ref.converged
    f_ref = 1.0 / ref.period

    got = {}
    for tones in ([1, 2, 3], [1, 3]):
        ps = ProbeShooting(fac(0.2), 'v', npts=250, tones=tones)
        f, amps, ph, info = ps.solve_multitone(
            f_lc, [2.0] + [0.05] * (len(tones) - 1), tol=1e-7, maxiter=12)
        got[tuple(tones)] = (f, info['converged'], abs(f - f_ref) / f_ref)

    full = got[(1, 2, 3)]
    pruned = got[(1, 3)]
    assert full[1] and pruned[1], \
        'both must CONVERGE -- the point is that convergence does not ' \
        'distinguish them: %r' % (got,)
    assert full[2] < pruned[2], \
        'including the even tone must be more accurate on an asymmetric ' \
        'circuit: full %.3e against pruned %.3e' % (full[2], pruned[2])
    assert pruned[2] / full[2] > 2.0, \
        'the pruning penalty should be substantial, got only %.2fx' \
        % (pruned[2] / full[2])


@pytest.mark.slow
def test_the_pac_probe_jacobian_agrees_with_finite_difference_and_is_cheaper():
    """`dI/dV` from K LINEAR PAC solves instead of 2K nonlinear ones.

    Measured on van der Pol, same fixture as the multitone solve::

        K  route   f          df/f        solves   wall
        2  FD      0.159124   +5.921e-02    16      14.2s
        2  PAC     0.159124   +5.921e-02     9      13.7s
        3  FD      0.150167   -4.099e-04    43      40.6s
        3  PAC     0.150167   -4.099e-04    19      29.8s

    Identical frequency to every printed digit, 2.3x fewer solves at K=3, and
    the gain GROWS with K because FD is O(2K) nonlinear solves against one
    nonlinear plus O(K) linear.

    ⚠⚠ FOUR DEFECTS STOOD BETWEEN "PAC HAS THE RIGHT QUANTITY" AND A WORKING
    JACOBIAN, and every one was found by a STRUCTURED discrepancy rather than
    by fitting a constant -- which is why none of them was papered over:

      * ratios of exactly 1, 2, 3 at K = 1, 2, 3  ->  `VS.vac` DEFAULTS TO 1,
        so every probe in the series chain was excited at once and, sharing one
        branch current, contributed K times over;
      * summing the folded pair CANCELLED and taking the larger HALVED  ->  the
        two entries at each harmonic SUBTRACT;
      * a sign-only error at m=3 with correct magnitude  ->  `PAC.solve` folds
        the sideband index away, so the pair's array ORDER is not stable across
        harmonics.  Ordering by magnitude passed at the solution (1.6e-05) and
        FAILED at the Newton's start (1.763): a heuristic that passes its gate
        and then fails in use.  Fixed by exciting at `j*f0 + delta`, which
        separates the pair in FREQUENCY -- direct at `m*f0 + delta`, image at
        `m*f0 - delta` -- so the rule is derived, not guessed;
      * the solve diverging to f = 0.0348 while `pac_jacobian` validated at
        1e-04  ->  a CHAIN RULE was missing.  The Jacobian is
        `d(Re I, Im I)/d(Re V, Im V)`; the unknowns are `(A, phi)` in degrees.
        A correct derivative wired to the wrong variables.

    ⚠ The last one is the reason `use_pac` re-validates at the Newton's OWN
    starting point rather than trusting a standalone check: a correct Jacobian
    and a broken solve coexisted, and only validating in situ caught it.
    """
    import warnings as _w
    from pycircuit.circuit.shooting import ProbeShooting
    circuit.default_toolkit = circuit.numeric

    def fac():
        def build():
            c = SubCircuit()
            c.add_node('v')
            c.add_node('x')
            c['C'] = C('v', gnd, c=1.0)
            c['RL'] = R('v', 'x', r=1e-2)
            c['L'] = L('x', gnd, L=1.0)
            c['B'] = BSource('v', gnd, gnd, 'v',
                             i_func=lambda u: (u - u ** 3 / 3.0))
            return c
        return build

    f_lc = 1.0 / (2.0 * np.pi)

    ## (1) It validates AT the solution and AWAY from it.  The second is the
    ## case the magnitude-ordering heuristic failed, so it is the one that
    ## matters -- a Jacobian is used away from the solution by definition.
    for amps, where in (([2.0, 0.0, -0.25], 'solution'),
                        ([2.0, 0.05, 0.05], 'Newton start')):
        ps = ProbeShooting(fac(), 'v', npts=300, harmonics=3)
        J, info = ps.pac_jacobian(f_lc, amps, [90.0] * 3, validate=True)
        assert info['validated'], where
        assert info['validation_reldiff'] < 1e-3, \
            'PAC vs FD at the %s: %.3e' % (where, info['validation_reldiff'])
        assert J.shape == (6, 6)

    ## (2) ⚠ THE OFF-DIAGONALS ARE THE ONLY ENTRIES THAT TEST THE SIDEBAND MAP.
    ## At K=1 there are none (m = j, so k = 0), which is exactly why a passing
    ## K=1 check missed the `vac` defect for three iterations of debugging.
    ps1 = ProbeShooting(fac(), 'v', npts=250, harmonics=1)
    J1, i1 = ps1.pac_jacobian(f_lc, [1.99], [90.0], validate=True)
    assert i1['validated'] and J1.shape == (2, 2)

    ## (3) END TO END: the PAC route must reach the SAME orbit as FD, since a
    ## faster wrong answer is worthless.
    got = {}
    for use_pac in (False, True):
        ps = ProbeShooting(fac(), 'v', npts=300, harmonics=3)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            f, amps, ph, info = ps.solve_multitone(
                f_lc, [2.0, 0.05, 0.05], tol=1e-7, maxiter=12, use_pac=use_pac)
        assert info['converged'], 'use_pac=%s did not converge' % use_pac
        got[use_pac] = (f, info['evaluations'])
    f_fd, n_fd = got[False]
    f_pac, n_pac = got[True]
    assert abs(f_pac - f_fd) / f_fd < 1e-6, \
        'the PAC route found a different orbit: %.8f against %.8f' % (f_pac, f_fd)
    assert n_pac < n_fd, \
        'the PAC route must cost fewer solves, got %d against %d' % (n_pac, n_fd)


def test_injection_locking_range_matches_averaging_on_the_textbook_control():
    """A6, first gate: the locking range of a driven van der Pol, read as the
    Floquet stability edge of the forced orbit under continuation, against
    the averaging prediction `dw_lock = I_inj / (2 C A)` -- NAMED before it
    was measured.

    Plain van der Pol, C = L = 1, a = 0, I_inj = 0.2 mu A.  Measured
    2026-09-07: the locked branch is tracked to exactly 1.00x the prediction
    (|lam| 0.98858 -> 0.99856, amplitude 2.176 -> 2.002 V) and lost at 1.10x.
    The saddle-node on the invariant circle, presenting as it should.

    ⚠⚠ THREE THINGS THIS TEST IS BUILT NOT TO DO, each of which produced a
    plausible wrong answer on the way here:
      * read lock from CONVERGENCE -- every detuning converges;
      * sweep from a FIXED seed -- it hops branches;
      * use a Q-based Adler formula -- three Q conventions for this circuit
        (32 / ~100 / ~201) give predictions 6x apart, and the measured edge
        sat near one of them by accident.  The averaging form has no Q.
    And my first injection scale used the TANK current (200x the
    negative-resistance current), a 10x overdrive that gave a smooth,
    monotone, wrong curve.

    This is the CONTROL: averaging is textbook here.  The hostile fixture
    (asymmetric, non-unit-reactance) is gated separately against the number
    this one validates.
    """
    last, pred_hz, f0, trace = _injection_lock_edge(_vdp_injected())
    assert last is not None, 'no locked branch found at all: %r' % (trace[:3],)
    assert 0.85 <= last <= 1.15, \
        'locked branch tracked to %.2fx the averaging prediction ' \
        '(I_inj/(2CA) = %.4e Hz, %.3f%% of f0); measured 1.00x. Trace: %r' % (
            last, pred_hz, 100 * pred_hz / f0,
            [(r, round(am, 3), round(lm, 5), lk) for r, am, lm, lk in trace])
    ## and the multiplier must have RISEN toward 1 along the branch -- the
    ## signature of the saddle-node, not of a branch that merely stopped
    lams = [lm for r, am, lm, lk in trace if lk]
    assert lams[-1] > lams[0] and lams[-1] > 0.995, \
        '|lam| along the locked branch went %.5f -> %.5f; it must rise ' \
        'toward 1 at the edge' % (lams[0], lams[-1])


def _pll_loop(K, kvco, fref=1e6, df=0.0):
    """VCO + multiplier PD + RC filter, `f0` EXACTLY `fref` and modulus 1 --
    so at zero gain the phase advances one cycle per reference period and
    EVERY offset is a solution, which is the marginal mode itself."""
    from pycircuit.circuit.elements_hdl import VcoHdl
    c = SubCircuit()
    for nd in ('ref', 'vco', 'ph', 'pd', 'ctl'):
        c.add_node(nd)
    c['Vref'] = VSin('ref', gnd, va=1.0, freq=fref, phase=0.0)
    c['X1'] = VcoHdl('ctl', gnd, 'vco', gnd, 'ph', f0=fref + df, kvco=kvco,
                     va=1.0, modulus=1.0)
    c['PD'] = _PllMultPd('vco', gnd, 'ref', gnd, 'pd', gnd, k=K)
    c['Rf'] = R('pd', 'ctl', r=1e3)
    c['Cf'] = C('ctl', gnd, c=1e-9)
    return c


def _pll_lambda(K, kvco, offset=0.0, fref=1e6, npts=400, df=0.0):
    """`|lambda|_max` of the closed loop, seeded at `offset` CYCLES on the
    integrator's own accumulator.  Returns None if the solve does not
    converge."""
    import warnings as _w
    cir = _pll_loop(K, kvco, fref, df)
    names = [str(n) for n in cir.nodes if str(n) != 'gnd!']
    x0 = np.zeros(cir.n - 1)
    ## ⚠ the ACCUMULATOR, not the `ph` node: `ph` is a dependent output and
    ## seeding it leaves the integrator's state untouched, so the solve
    ## returns to the same equilibrium (measured: offsets on `ph` alone
    ## never left the saddle).
    x0[names.index('X1._state0')] = offset
    x0[names.index('ph')] = offset
    pss = PSS(cir, method='gear', reltol=1e-10)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=1.0 / fref, timestep=1.0 / fref / npts, x0=x0,
                  maxiterations=100)
    if not pss.converged:
        return None
    mu = np.linalg.eigvals(np.asarray(pss._monodromy))
    return float(np.max(np.abs(mu)))


def test_closing_a_pll_loop_pins_the_marginal_phase_mode_at_the_loop_bandwidth():
    """A6 step 2, and it decides whether the rest of A6 is a build or research.

    The roadmap left this open: "a free-running integrator's phase row is a
    marginal mode -- measured dx_end/dx_0 = 1.000000 -- so `I - M` is singular
    there by construction.  Locking it needs feedback (a PLL) ... the fold
    repairs the residual's VALUE, not the Jacobian's RANK."  So: does closing
    the loop move that multiplier off 1?

    IT DOES, at exactly the loop bandwidth.  Linearising the averaged loop --
    equilibrium at `f0 + kvco Vctl = fref`, `Vctl = (K/2) cos(2 pi theta)`,
    slope `-+ pi K` V/cycle -- gives `dtheta/dt = -+ kvco pi K theta`, so over
    one reference period

        |lambda|_phase = exp( -+ pi * kvco * K * T )

    MEASURED 2026-09-18 to 4e-09..4e-05 relative over two decades of gain, on
    TWO independent knobs (K at fixed kvco, and kvco at fixed K), which agree
    to all printed digits wherever the PRODUCT `kvco*K` matches -- so it is the
    loop-gain product that governs and nothing else.

    ⚠ BOTH EQUILIBRIA ARE REAL AND THE SOLVER DOES NOT PREFER THE STABLE ONE.
    A multiplier PD has two per cycle, `theta = +-1/4`: measured 0.999371484
    (stable) and 1.000628121 (saddle) at K = 1e-2, kvco = 2e4.  From the
    NATURAL zero seed the solve lands on the SADDLE -- so convergence does not
    select stability, and lock must be read from the multiplier.  That is the
    same lesson A6's injection-locking gates already carry ("read lock from
    CONVERGENCE -- every detuning converges").

    ⚠ The branch is selected by the integrator's ACCUMULATOR, not by the `ph`
    node, which is a dependent output.
    """
    fref, T = 1e6, 1e-6
    kvco, K = 2e4, 1e-2
    want_stable = np.exp(-np.pi * kvco * K * T)
    want_saddle = np.exp(+np.pi * kvco * K * T)

    ## the NEGATIVE CONTROL: with the loop open the driven solve is
    ## underdetermined -- every phase offset is a solution, the solve returns
    ## the seed's, and its phase multiplier sits at 1 (measured 1 + 1.8e-11).
    ## ⚠ THIS USED TO ASSERT NON-CONVERGENCE, and that held only because the
    ## seed (offset 0) put the orbit's start on the phase WRAP: the Newton
    ## iterates straddled it and stalled.  Off the wrap the open loop
    ## converged on that code too (offset 0.3: 1.000000000018), and since
    ## `_wrap_jump` (2026-09-23) it converges on it as well.  Singular
    ## `I - M` is a family of orbits, not a failed solve; the multiplier
    ## says which.
    lam_open = _pll_lambda(0.0, kvco, offset=0.3)
    assert lam_open is not None and abs(lam_open - 1.0) < 1e-9, \
        'with K = 0 the phase mode must be marginal: |lambda| = %r' % lam_open

    lam_saddle = _pll_lambda(K, kvco, offset=0.0)
    lam_stable = _pll_lambda(K, kvco, offset=0.25)
    assert lam_saddle is not None and lam_stable is not None
    assert abs(lam_saddle / want_saddle - 1.0) < 1e-5, (lam_saddle, want_saddle)
    assert abs(lam_stable / want_stable - 1.0) < 1e-5, (lam_stable, want_stable)
    ## and they must straddle 1 -- feedback PINS the mode, it does not merely
    ## perturb it in some direction
    assert lam_stable < 1.0 < lam_saddle, (lam_stable, lam_saddle)

    ## the PRODUCT is what governs: ten times the gain on the other knob is
    ## the same loop, to the digit.
    a = _pll_lambda(K, kvco, offset=0.0)
    b = _pll_lambda(K * 10.0, kvco / 10.0, offset=0.0)
    assert abs(a / b - 1.0) < 1e-9, \
        'kvco*K is the loop gain, so these must agree: %.12f vs %.12f' % (a, b)


def test_the_pll_lock_range_is_the_loop_bandwidth_and_the_whole_locus_is_predicted():
    """A6 step 3: the hold-in range, gated as a LOCUS rather than an edge, and
    with every constant imported from the two measurements before it.

    For a first-order loop the VCO can be pulled by at most `kvco*(K/2)`, so the
    hold-in range is `|f0 - fref| <= Df_max = kvco*K/2` -- THE SAME CONSTANT
    step 2 read off the Floquet multiplier and step 4 read out of the noise
    corner.  Inside it the equilibrium sits at `cos(2 pi theta) = -Df/Df_max`,
    the linearised rate is `kvco pi K |sin(2 pi theta)|`, and with step 4's
    second pole:

        -ln|lam| / (pi kvco K T) = s + s^2 (f_c/f_RC),   s = sqrt(1 - x^2)

    a unit semicircle in `x = Df/Df_max`, tilted by the RC.  NOTHING IS FITTED:
    `Df_max` and `f_c` come from `kvco` and `K`, and `f_RC` from R and C.
    MEASURED 2026-09-18 to 6.3e-07 .. 7.9e-06 over x = 0 .. 0.99, with lock lost
    between 99 and 100.5 Hz.  So this is a THIRD independent route to one
    constant -- a multiplier, a noise corner, and a detuning edge.

    ⚠ Seeded on the STABLE branch and continued, never from a fixed seed: step 2
    showed the natural seed lands on the saddle, and A6's own injection-locking
    record says a fixed seed hops branches.
    """
    fref, T = 1e6, 1e-6
    kvco, K, rf, cf = 2e4, 1e-2, 1e3, 1e-9
    dfmax = kvco * K / 2.0
    f_rc = 1.0 / (2.0 * np.pi * rf * cf)

    for df in (0.0, 60.0, 90.0, 99.0):
        lam = _pll_lambda(K, kvco, offset=0.25, df=df)
        assert lam is not None and lam < 1.0, (df, lam)
        x = df / dfmax
        s = np.sqrt(1.0 - x ** 2)
        want = s + s * s * (dfmax / f_rc)
        got = -np.log(lam) / (np.pi * kvco * K * T)
        assert abs(got / want - 1.0) < 1e-4, (df, got, want)

    ## the EDGE: beyond the hold-in range there is no locked branch at all
    assert _pll_lambda(K, kvco, offset=0.25, df=105.0) is None, \
        'a locked branch survived past Df_max = %.1f Hz' % dfmax
