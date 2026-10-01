"""Shooting tests: shooting ppv.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.circuit.simwarnings import AccuracyWarning, UsageWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
import warnings
from pycircuit.circuit.tests._shooting_fixtures import (_adjoint_ladder,
    _comparator_relaxation_oscillator,
    _exact_relaxation_oscillator_ppv,
    _injection_lock_edge,
    _loss_osc,
    _orbit_modulated_vdp,
    _raw_pair_integrals,
    _relaxation_oscillator_seed,
    _solve_slow,
    _vdp_asym,
    _vdp_injected,
    _vdp_with_noise)


def _vdp_ppv(npts, mu=1.0):
    """A converged van der Pol and its PPV."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=6.6634, timestep=6.6634 / npts,
                  x0=np.array([2.0, 0.0]), maxiterations=50)
    assert pss.converged, 'van der Pol did not converge at %d points' % npts
    v, info = pss.ppv()
    return c, pss, v, info


def test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has():
    """⚠ THE GATE FOR A2, and it is a physical experiment, not an identity.

    A PPV is only worth having if `v . delta` is the phase shift a real
    perturbation causes. So: displace van der Pol's state by `eps delta`,
    integrate the FULL NONLINEAR system for ONE period, and read the
    displacement along the orbit tangent. Nothing in that measurement
    touches the monodromy, the adjoint, or the bordered solve.

    ⚠⚠ IT USED TO INTEGRATE THREE PERIODS "UNTIL THE TRANSVERSE
    COMPONENTS HAVE DIED", AND THAT REASON WAS WRONG -- MEASURED. The
    premise is false at high Q (`lambda_2^3 = 0.954` at Q = 64) and the
    waiting was not doing the work at ANY Q:

        npts   nper=1     nper=2     nper=3
         200   2.305e-02  2.348e-02  2.346e-02
         400   1.058e-02  1.106e-02  1.104e-02
         800   4.889e-03  5.370e-03  5.377e-03

    ⚠ `nper` BARELY MATTERS AND ONE PERIOD IS SLIGHTLY BETTER, while the
    error halves per `npts` doubling.  The residual is O(h)
    DISCRETISATION of the map, not transverse contamination -- so there
    was nothing to wait for, and the wait cost 3x for nothing.

    ⚠ CONFIRMED FROM THE OTHER END: at `mu = 1` the contamination is
    EXACTLY zero (`lambda_2^3 = 6e-10`) and the error is still 6.8% at a
    transverse ratio of 3.  A 6.8% error with zero contamination is
    larger than the contamination itself at Q = 64 (2.3%).  Contamination
    was never the dominant term in this gate.

    ⚠⚠ AND THREE EXTRAPOLATION REPAIRS WERE MEASURED AND ALL MADE IT
    WORSE, recorded so nobody re-derives them:

        Q      raw n=1   raw n=3   Aitken   lambda_2-extrapolation
        0.12   6.761%    6.898%    7.038%   6.901%
        15.9   1.301%    1.623%    3.259%   4.074%
        63.7   2.314%    2.520%    5.760%   9.039%

    Aitken and the `lambda_2` extrapolation both assume the residual IS
    the geometric transverse mode.  It is not, so they amplify what is
    left instead of removing it -- and the `lambda_2` form amplifies by
    `1/(1 - lambda_2)`, which is `Q`, so it degrades fastest exactly
    where it was supposed to help.

    ⚠⚠ THAT TABLE IS ONE SIDE OF A CROSSOVER, NOT THE GENERAL CASE
    (the docs session's criterion, 2026-09-09, MEASURED here at m = 20).
    Aitken helps iff the geometric transverse part of the raw error
    exceeds the O(h) floor; the partition is read BEFORE Aitken enters,
    by whether the raw n=3 error moves with `npts`.  This gate sits on
    the floor-dominated side (a random 2-D kick is ~71 % tangential, the
    raw error halves per doubling), so Aitken loses.  On the same van
    der Pol with 18 RC branches (m = 20, one slow branch at 100 T, kicks
    ~0.2 tangential, scratchpad aitken_m20.py):

        kick        raw n=3 at npts 100/200/400   Aitken(1,2,3)   scored
        random0     0.499  0.459  0.456  (flat)   0.178  0.052    Aitken wins 2.8x / 8.8x
        random1     0.735  0.695  0.681  (flat)     --   0.055    Aitken wins 12.7x
        transverse  1.56   1.64   1.58   (flat)   1.34   0.084    Aitken wins
        random2     0.071  0.038  0.057  (moves)  0.134  0.049    raw wins
        tangent     8.8e-3 2.8e-3 7.5e-4 (halves) 0.18   0.12     raw wins (control)

    Flat rows: Aitken wins; the moving row and the tangent control: raw
    wins.  Rows whose successive-difference ratios differ by more than
    10 % (18 packed slow modes) are NOT a geometric sequence and test
    nothing; at npts = 400 every row read that way until the transient's
    `reltol` went from 1e-9 to 1e-11 -- the adaptive integrator's own
    endpoint error, ~5e-3 of the signal, once the geometric differences
    had shrunk to meet it.  ⚠ A constant-sequence control (the tangent)
    cannot see an n-dependent error, so it does not certify the others.
    Two floors, two knobs (the docs session's refinement of its own
    criterion): the PSS grid's O(h), which refining `npts` LOWERS, and
    the transient's `reltol`, which refining `npts` UNCOVERS.  Before
    scoring a row non-geometric, `reltol` must sit ~2 decades below the
    SMALLEST difference extrapolated (in state units, `eps` times the
    difference), not below the signal: Aitken differences twice and eats
    two decades of headroom.  Here the differences were ~1e-8 of state
    against 2e-9 bought by reltol 1e-9 -- under one decade.  Aitken removes ONE
    mode, so its gain is bounded by the spread of the transverse rates,
    not by contamination/floor -- 3x to 20x here, 350x at m = 2 where
    the leftover is the floor.  So "the repairs made it worse" is a fact
    about THIS gate, and a real circuit's random kick (mostly transverse)
    is on the other side of the crossover.

    ⚠ THE SCALE IS THE ASSERTION, NOT JUST THE DIRECTION. A direction check
    passes for any normalisation, and the normalisation is exactly what was
    in doubt. Measured against the true shift, per doubling of the period
    grid:

        npts   worst |1-ratio|   rel resid   fitted scale
         200      5.20e-02        2.305e-02    1.003584
         400      2.31e-02        1.058e-02    1.001969
         800      1.05e-02        4.889e-03    1.000978

    The fitted scale converges to ONE -- not to some other constant that a
    direction-only test would have accepted -- and the residual falls at
    O(h), which is the discretisation of the map the PPV is computed from.

    ⚠ AND THE PREMISE ABOUT DECAY IS FALSE AT HIGH Q, WHICH THIS GATE
    SURVIVES FOR A REASON THAT IS NOT ITS OWN.  Three periods kills the
    transverse mode at `mu = 1` (`lambda_2^3 = 6e-10`) and does nothing at
    `Q = 64` (`lambda_2^3 = 0.954`).  MEASURED there anyway: the raw error
    stays under 1% -- 0.73% at n=1 rising to 0.88% at n=6 -- because a
    random direction in TWO dimensions is ~71% tangential, so the phase
    signal dominates and the contamination is a small additive term.

    ⚠⚠ THAT PROTECTION SCALES AS `1/sqrt(m)` AND VANISHES ON A REAL
    CIRCUIT.  The tangential fraction of a random unit vector is
    `~1/sqrt(m)`: 0.71 at m=2, 0.32 at m=10, 0.10 at m=100.  So on any
    circuit with more than a handful of states a random kick is MOSTLY
    TRANSVERSE and the three-period premise does bite.  **This gate is
    sound at m = 2 and would not be at m = 20, with nothing in it
    changing.**  The transverse ratios asserted below are what makes that
    visible instead of implicit.
    """
    from pycircuit.circuit.transient import Transient
    rng = np.random.default_rng(0)
    dirs = [d / np.linalg.norm(d) for d in
            (rng.standard_normal(2) for _ in range(4))]

    out = []
    ratios = []
    for npts in (200, 400):
        cir, pss, v, info = _vdp_ppv(npts)
        m = cir.n - 1
        irn = pss.irefnode
        T = pss.period
        xdot = info['xdot']
        x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))

        def integrate(xi, nper=1, ppp=2000):
            ## ⚠ ONE PERIOD, NOT THREE -- see the docstring. Waiting for
            ## the transverse mode was never what made this work, and
            ## dropping it is 3x cheaper AND marginally more accurate.
            ## ⚠ `reltol` HERE IS THE COST, not `ppp`. `Transient` adapts,
            ## so `timestep` is a first step and the tolerance sets the
            ## step count. At 1e-12 this test took 447 s; the signal being
            ## measured is ~9e-3, so 1e-9 is still three orders finer than
            ## anything it has to resolve.
            tran = Transient(cir, toolkit=circuit.numeric, reltol=1e-9,
                             iabstol=1e-13, vabstol=1e-11)
            with quiet():
                res = tran.solve(refnode=gnd, tend=nper * T,
                                 timestep=T / ppp, x0=xi)
            return np.asarray(res.x, dtype=float)[:, -1]

        ref = integrate(x0f)
        meas, pred = [], []
        for d in dirs:
            eps = 1e-5
            dr = np.concatenate((d[:irn], np.zeros(1), d[irn:]))
            dx = np.delete(integrate(x0f + eps * dr) - ref, irn)
            meas.append(float(dx @ xdot) / float(xdot @ xdot) / eps)
            pred.append(float(v[:m] @ d))
        ## ⚠ HOW TRANSVERSE EACH KICK ACTUALLY IS, asserted rather than
        ## assumed.  Signal and contamination are the SAME KNOB: a
        ## tangential displacement IS a pure phase shift, so it has
        ## nothing transverse to decay, and the direction that maximises
        ## `v.d` is exactly the one that tests transverse decay LEAST.
        ## Measured over all directions at Q = 15.9: `|v.d|` spans 1.6e-4
        ## to 5.0e-1, and the transverse ratio runs the opposite way --
        ##
        ##     |b/a|      signal      contamination at n = 3
        ##      0          max         0
        ##      1          1.4x down   1.0%
        ##      3          3.4x down   3.4%
        ##     40         40x down    41%
        ##
        ## ⚠ SO A GATE CANNOT MAXIMISE BOTH, and one that reports a
        ## beautiful agreement may simply have kicked along the orbit.
        ## These four seeded directions give |b/a| = 0.58, 14.5, 0.92,
        ## 2.39 at mu = 1 and 0.92, 6.81, 1.43, 1.43 at Q = 15.9 -- the
        ## useful middle, by luck of the seed until this assertion.
        xh = xdot / np.linalg.norm(xdot)
        for d in dirs:
            tan = abs(float(d @ xh))
            ratios.append(np.sqrt(max(1.0 - tan ** 2, 0.0)) / max(tan, 1e-300))
        meas, pred = np.array(meas), np.array(pred)
        out.append((np.linalg.norm(meas - pred) / np.linalg.norm(meas),
                    float((pred @ meas) / (pred @ pred))))

    (r_coarse, _s_coarse), (r_fine, s_fine) = out
    assert abs(s_fine - 1.0) < 0.02, \
        'the PPV predicts the phase shift up to a factor of %.4f, not 1. ' \
        'The direction is right and the NORMALISATION is not -- which is ' \
        'the whole question this test exists to settle' % s_fine
    ratio = r_coarse / r_fine
    assert 1.6 < ratio < 2.8, \
        'the disagreement with the measured phase shift falls %.2fx per ' \
        'doubling (%.3e -> %.3e), not the O(h) that says it is ' \
        'discretisation of the period map' % (ratio, r_coarse, r_fine)
    ## ⚠ THE GATE NOW STATES HOW MUCH TRANSVERSE DECAY IT EXERCISES.
    ## Without this it could pass having kicked almost along the orbit,
    ## where there is nothing to decay -- which is exactly how a
    ## well-conditioned-looking probe tests nothing.
    ratios = np.array(ratios)
    assert ratios.max() > 1.0, \
        'every direction is more tangential than transverse (max |b/a| = ' \
        '%.3f); this gate would pass without exercising transverse decay ' \
        'at all' % ratios.max()
    assert ratios.min() < 3.0, \
        'every direction is strongly transverse (min |b/a| = %.3f), so ' \
        'the signal is far down and the comparison is a difference of ' \
        'near-zeros' % ratios.min()

def test_the_ppv_normalisation_is_not_the_one_transcribed():
    """⚠ `v . q = 1` IS THE WRONG SCALE HERE, and it looks right.

    Demir's Remark 3.1 reads `v_1^T C(0) u_1(0) = 1`, and bordering the
    augmented system with `q = C xdot` makes `v . q = 1` fall out for free.
    On this fixture `v . q = -1.0696` where `v . xdot = 1`: a factor
    2.07 and the opposite sign (an earlier version of this docstring
    said "7%"; that was `|v . q| - 1`, not the error -- corrected by the
    review session's audit).  The vector this bordered solve returns behaves as
    `C^T v_1` -- it contracts with a STATE perturbation directly -- so the
    two statements are both true of different objects, and using one where
    the other belongs is a silent scale error in every phase-noise number
    downstream.

    ⚠ THE OBVIOUS REPAIR IS ALSO WRONG, and was measured so before this
    normalisation was chosen: predicting the shift as `v^T C delta`
    gives residuals of 0.36 / 0.40 / 0.42 that GROW under refinement, with
    per-direction ratios scattering from -0.44 to 28.7. `v . delta` with
    `v . xdot = 1` converges at O(h). Do not "fix" this back without
    re-running that experiment.
    """
    cir, _pss, v, info = _vdp_ppv(400)
    m = cir.n - 1
    assert abs(float(v[:m] @ info['xdot']) - 1.0) < 1e-9, \
        'the returned PPV is not normalised by v . xdot = 1'
    vq = float(v[:m] @ info['q'])
    assert abs(vq - 1.0) > 0.5, \
        'v . q is %.6f, i.e. indistinguishable from the normalisation this ' \
        'test exists to rule out. Either C is now the identity on this ' \
        'circuit (in which case pick another) or the scale has been ' \
        'changed back' % vq


def test_the_ppv_border_residual_is_the_free_check():
    """`y` coming back zero is D&R's own correctness check, and it is real.

    With a zero first block on the right-hand side, `(I - M^T)v + y q = 0`
    forces `y q = 0`, so a nonzero `y` means the border absorbed a residual
    that belongs to the null space -- the computed `v` is not in it. Both
    solves report it, and both come back at ~1e-11.
    """
    _cir, _pss, _v, info = _vdp_ppv(200)
    for key in ('border_residual', 'tangent_border_residual'):
        assert abs(info[key]) < 1e-7, \
            '%s is %.3e; the bordered system absorbed a null-space ' \
            'residual, so the vector it returned is not the PPV' \
            % (key, info[key])
    assert info['null_residual'] < 1e-7, \
        '||v - M^T v|| / ||v|| is %.3e' % info['null_residual']


def test_the_ppv_refuses_what_has_no_phase():
    """A driven circuit's phase is its source's, not its own."""
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    assert pss.converged and not pss.autonomous
    with pytest.raises(ValueError, match='FREE-RUNNING'):
        pss.ppv()


def test_the_ppv_samples_are_the_propagated_adjoint():
    """`v(t)` over the whole period, checked by FORWARD machinery.

    Oscillator phase noise is an integral over the orbit — Demir's
    diffusion constant is `c = (1/T) ∫ v₁ᵀ(t) B(t)Bᵀ(t) v₁(t) dt` — so the
    PPV at `t = 0` is not enough. `Phi(T,s)ᵀ v(T) = v(s)`, and the reverse
    replay computes that sequence on its way to the answer; it was being
    discarded.

    ⚠ CHECKED AGAINST THE FORWARD STEP MAPS, not against itself. The
    propagator `Phi(T,s_j)` is rebuilt here by applying the forward
    recursion from step `j` to the end, and its transpose is compared with
    what the reverse pass recorded. Same shape as the existing `M` vs `Mᵀ`
    test, but per step, and it is what a future oscillator pnoise would be
    resting on. Measured worst case over all steps: 1.8e-15.
    """
    _cir, pss, v, info = _vdp_ppv(60)
    fp = pss.factored_period()
    n = fp.width
    m = pss.cir.n - 1
    ## ⚠ the RAW pair adjoint: `samples` is now the pair-CONSISTENT
    ## contraction (see `ppv`), which is not `Psi^T v` and must not be.
    samples = np.asarray(info['samples_pair'])
    assert samples.shape == (len(fp.steps), n), \
        'expected one PPV sample per step, got %r' % (samples.shape,)

    cs0, cs1, ring = [], [], list(fp.opening)
    for _lu, C_new, _a, _b in fp.steps:
        cs0.append(ring[0])
        cs1.append(ring[1])
        ring = [C_new, ring[0]]

    def fwd_from(j, p):
        p0, p1 = p[:m].copy(), p[m:].copy()
        for k in range(j, len(fp.steps)):
            lu, _Cn, alphas, b = fp.steps[k]
            assert not b, 'this rebuild assumes the Gear-2 companion'
            pn = -lu.solve(alphas[1] * (cs0[k] @ p0)
                           + alphas[2] * (cs1[k] @ p1))
            p0, p1 = pn, p0
        return np.concatenate((p0, p1))

    ## the rebuild must be the monodromy, or it is checking nothing
    Mf = np.column_stack([fwd_from(0, e) for e in np.eye(n)])
    Mm = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    assert np.linalg.norm(Mf - Mm) / np.linalg.norm(Mm) < 1e-12, \
        'the forward rebuild is not the monodromy, so it cannot check the ' \
        'reverse pass'

    worst = 0.0
    for j in range(len(fp.steps)):
        Psi = np.column_stack([fwd_from(j, e) for e in np.eye(n)])
        pred = Psi.T @ v
        worst = max(worst, np.linalg.norm(pred - samples[j])
                    / max(np.linalg.norm(pred), 1e-300))
    assert worst < 1e-11, \
        'the reverse pass states are not the propagated adjoint (worst ' \
        '%.3e), so they are not the PPV over the period' % worst


def test_no_periodic_covariance_exists_for_an_oscillator():
    """⚠ WHY OSCILLATOR NOISE NEEDS THE PPV AND NOT A COVARIANCE SHOOT.

    Time-varying noise statistics are a Lyapunov ODE alongside the
    transient, and for a periodic large signal its periodic solution is a
    shooting problem whose monodromy is the KRONECKER SQUARE of the
    circuit's: `Phi_lyap = M ⊗ M`, so its multipliers are the pairwise
    products `lambda_i lambda_j`.

    For a DRIVEN circuit that is a single linear solve — the Lyapunov
    equation is linear in `K`, so no Newton iteration. For an AUTONOMOUS
    one it does not exist: `lambda_1 = 1` gives `lambda_1^2 = 1`, and
    `I - M ⊗ M` is exactly as singular as `I - M`. The covariance does not
    settle, it GROWS — variance linear in `t` is a random walk, which is
    phase diffusion, which is the linewidth.

    ⚠ SO THE UNIT MULTIPLIER IS NOT AN INCONVENIENCE HERE, IT IS THE
    ANSWER. The same near-unit-eigenvalue obstruction this codebase keeps
    meeting — eigen-selection, the bordered phase row, the PPV — appears
    once more as `lambda_1^2 = 1`, and this time what it obstructs is the
    wrong method for the question. Relayed measurements on three LTP
    systems put the kron identity at 2.2e-14 and the phase-mode growth at
    dead linear (trace K 1.825 -> 1793.4 over 1000 periods); this checks
    the structural half on our own monodromy.
    """
    circuit.default_toolkit = circuit.numeric

    _cir, pss, _v, _info = _vdp_ppv(60)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    K = np.kron(M, M)

    ev = np.linalg.eigvals(M)
    outer = np.sort_complex(np.array([a * b for a in ev for b in ev]))
    got = np.sort_complex(np.linalg.eigvals(K))
    scale = max(float(np.max(np.abs(outer))), 1e-30)
    assert np.linalg.norm(got - outer) / scale < 1e-9, \
        'the covariance monodromy is not the Kronecker square of the ' \
        'circuit monodromy, so its multipliers are not the pairwise products'

    s_lyap = np.linalg.svd(np.eye(n * n) - K, compute_uv=False)[-1]
    s_circ = np.linalg.svd(np.eye(n) - M, compute_uv=False)[-1]
    assert s_lyap < 1e-7, \
        'I - M kron M has sigma_min %.3e on an AUTONOMOUS circuit. If the ' \
        'covariance shooting problem has become solvable, the unit ' \
        'multiplier has gone, and so has the oscillator' % s_lyap
    assert s_lyap < 100 * max(s_circ, 1e-16), \
        'the covariance system is singular to a different degree than the ' \
        'circuit one (%.3e against %.3e); the obstruction should be the ' \
        'same unit multiplier squared' % (s_lyap, s_circ)

    ## the contrast: a DRIVEN circuit has no such obstruction
    cir2 = _adjoint_ladder(4)
    pss2 = PSS(cir2, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    fp2 = pss2.factored_period()
    n2 = fp2.width
    M2 = np.column_stack([fp2.matvec(e) for e in np.eye(n2)])
    s2 = np.linalg.svd(np.eye(n2 * n2) - np.kron(M2, M2),
                       compute_uv=False)[-1]
    assert s2 > 1e-3, \
        'I - M kron M is near-singular (%.3e) on a DRIVEN circuit too, so ' \
        'the obstruction is not the unit multiplier after all' % s2


@pytest.mark.parametrize('mu', [0.5, 1.0])
def test_the_coloured_noise_functional_is_exactly_zero_here(mu):
    """⚠ A REGRESSION TEST WITH AN EXACT ANSWER OF ZERO, and it separates
    two functionals that share a PPV.

    Demir 2002 defines two different scalars from the same `v_1(t)`:

        c_w  = (1/T) ∫ v_1ᵀ B_w B_wᵀ v_1 dt     WHITE   — QUADRATIC
        V_0m = (1/T) ∫ v_1ᵀ B_cm dt             COLOURED — LINEAR

    The white one is the time-average of a SQUARE; the coloured one is the
    plain time-average — the zeroth Fourier coefficient of a periodic
    scalar. Using the quadratic form for a coloured source returns a
    plausible non-zero number from the same PPV, and nothing that did not
    know to look would catch it.

    §VIII gives the discriminating case: on a parallel-RLC oscillator with
    a nonlinear current source, "the time-average of [the Floquet vector
    entry] for the capacitor voltage is 0! Thus … any stationary …
    colored-noise source … connected across the capacitor has NO
    contribution to the oscillator spectrum due to phase noise, because
    V_0m = 0." Van der Pol is that circuit.

    ⚠ THE SECOND ASSERTION IS WHAT MAKES THIS A TEST RATHER THAN A
    TAUTOLOGY. A vector of zeros would pass the first one. The RMS of the
    same entries is ~0.40, so the WHITE functional is emphatically not
    zero in the same position — measured `|mean|/rms ~ 1e-11`. Zero and
    non-zero from one PPV, which is exactly the discrimination the
    coloured functional needs and the quadratic one destroys.

    (Both entries come back zero here, not just the capacitor's, which van
    der Pol's symmetry under `(v,i) -> (-v,-i)` predicts: the PPV is odd
    over the period.)
    """
    period = 6.6634 if mu >= 1.0 else 6.35
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=period, timestep=period / 200,
                  x0=np.array([2.0, 0.0]), maxiterations=60)
    assert pss.converged
    m = c.n - 1
    _v, info = pss.ppv()

    S = np.asarray(info['samples'])[:, :m]
    tms = np.asarray(info['times'], dtype=float)
    h = np.diff(tms)
    T = pss.period
    mean = (S * h[:, None]).sum(axis=0) / T           # the COLOURED scalar
    rms = np.sqrt(((S ** 2) * h[:, None]).sum(axis=0) / T)   # ~ the WHITE one

    assert (rms > 0.1).all(), \
        'the PPV samples are ~zero (rms %s), so the zero below would be ' \
        'vacuous' % np.array2string(rms, precision=4)
    ratio = np.abs(mean) / rms
    assert (ratio < 1e-8).all(), \
        'the time-average of the PPV is %s (relative %s), not zero. Demir ' \
        'section VIII says a coloured source on this oscillator ' \
        'contributes NO phase noise; a non-zero mean here means either ' \
        'the PPV or the averaging is wrong' \
        % (np.array2string(mean, precision=4),
           np.array2string(ratio, precision=4))


def test_the_null_residual_fires_on_a_wrong_ppv_and_says_how_blind_it_is():
    """A MUTATION check on `null_residual`, plus the bound it actually gives.

    Two sessions have now read the SAME flat `1e-9`/`1e-11` residual in
    opposite directions -- one as proof the bordered PPV solve stays accurate
    as `lambda_2 -> 1`, this file as the instrument saying nothing.  Neither
    is right, and the difference is measurable, so it is measured here rather
    than argued.

    `null_residual` is `||v - M^T v|| / ||v||`.  Inject a 1% error into the
    converged `v`:

      * in a RANDOM direction it reads 1.65e-02 against a converged floor of
        4.6e-11 -- nine orders.  The gate is real and every assertion on it in
        this file can fail.  That is the half the roadmap had too pessimistic.
      * along the `lambda_2` LEFT-EIGENDIRECTION it reads `0.01*(1 - lam2)`
        exactly: 1.003e-02, 9.950e-05, 1.000e-06, 1.000e-08 at
        `lam2 = 0.000856, 0.990049, 0.999900, 0.999999`.  So the error the
        residual cannot exclude is `null_residual / (1 - lam2)`, and the
        blindness grows without bound as the circuit gets better.  That is
        the half a "residual remains at 1e-9" claim gets wrong.

    `info['null_residual_amplification']` ships the second factor so a caller
    can convert one number into the other.  Read together or neither.
    """

    floor, injected = 4.6e-11, 0.01
    for tau, lam2_want in ((None, 0.000856), (1e2, 0.990049), (1e4, 0.999900)):
        cir, pss = _solve_slow(tau)
        with quiet(AccuracyWarning):
            v, info = pss.ppv()
        v = np.asarray(v, dtype=float).ravel()
        fp = pss.factored_period()
        n = fp.width
        v = v[:n]

        ## The left eigenvectors of `M` are the eigenvectors of `M^T`; the
        ## `lam2` one is the direction the residual is least able to see.
        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float).ravel()
                             for e in np.eye(n)])
        lams, W = np.linalg.eig(M.T)
        keep = [i for i in range(n) if abs(lams[i] - 1.0) > 1e-6]
        j = max(keep, key=lambda i: np.real(lams[i]))
        lam2 = float(np.real(lams[j]))
        w2 = np.real(W[:, j])
        w2 = w2 / np.linalg.norm(w2)

        def resid(vv):
            mv = np.asarray(fp.matvec_transposed(vv), dtype=float).ravel()
            return (np.linalg.norm(vv - mv)
                    / max(float(np.linalg.norm(vv)), 1e-300))

        assert abs(lam2 - lam2_want) < 1e-5, \
            'tau/T=%r: lam2 is %.6f, fixture expects %.6f' % (
                tau, lam2, lam2_want)

        nv = float(np.linalg.norm(v))
        rng = np.random.default_rng(7)
        d = rng.standard_normal(n)
        d = d / np.linalg.norm(d)

        r_clean = resid(v)
        r_rand = resid(v + injected * nv * d)
        r_lam2 = resid(v + injected * nv * w2)

        ## 1. The gate FIRES: a generic error is caught far above the floor.
        assert r_clean < floor * 10, \
            'tau/T=%r: converged residual %.3e is above the floor' % (
                tau, r_clean)
        assert r_rand > 1e4 * r_clean, \
            'tau/T=%r: a %g random error moved the residual only %.3e -> ' \
            '%.3e; the assertions on this key would be decorative' % (
                tau, injected, r_clean, r_rand)

        ## 2. And it is BLIND by exactly `1 - lam2` in the worst direction.
        assert abs(r_lam2 / (injected * (1.0 - lam2)) - 1.0) < 0.02, \
            'tau/T=%r: residual along the lam2 direction is %.3e, the ' \
            '(1-lam2) scaling predicts %.3e' % (
                tau, r_lam2, injected * (1.0 - lam2))

        ## 3. The shipped amplification is that factor, so
        ##    `null_residual * amplification` is the error it cannot exclude.
        amp = info['null_residual_amplification']
        assert abs(amp * (1.0 - lam2) - 1.0) < 1e-6, \
            'tau/T=%r: amplification %.6e does not match 1/(1-lam2) %.6e' % (
                tau, amp, 1.0 / (1.0 - lam2))
        assert info['null_residual'] * amp < 1e-4, \
            'tau/T=%r: the residual admits a relative error of %.3e in v' % (
                tau, info['null_residual'] * amp)


def test_a_slow_node_degrades_the_ppv_border_silently():
    """⚠ THE BORDER FIXES THE PHASE MODE AND NOTHING ELSE.

    `ppv()`'s bordered solve removes the singularity the UNIT multiplier
    causes. A second multiplier approaching 1 — which one slow node puts
    there — is untouched, and the conditioning goes with it. Measured on
    van der Pol with one weakly coupled RC node:

        tau/T   |lambda_2|   sigma_min(bordered)   null residual
        none     0.000856         8.62e-01            4.1e-11
        1e2      0.990049         4.49e-03            4.6e-11
        1e4      0.999900         4.47e-05            4.6e-11
        1e6      0.999999         4.47e-07            4.4e-11

    ⚠ `sigma_min` tracks `T/tau` over six decades WHILE THE RESIDUAL DOES
    NOT MOVE. GMRES converges, every diagnostic reads clean, and six
    digits of conditioning are gone. That is why `ppv()` estimates
    `|lambda_2|` explicitly instead of trusting a small residual — a
    residual cannot see this, and neither could this test if it only
    checked residuals.
    """
    smins, lam2s = [], []
    for tt in (1e2, 1e4):
        cir, pss = _solve_slow(tt)
        m = cir.n - 1
        fp = pss.factored_period()
        n = fp.width
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        with quiet(AccuracyWarning):
            _v, info = pss.ppv()
        q = info['q']
        qp = np.concatenate((q, np.zeros(n - m)))
        A = np.zeros((n + 1, n + 1))
        A[:n, :n] = np.eye(n) - M.T
        A[:n, n] = qp
        A[n, :n] = qp
        smins.append(np.linalg.svd(A, compute_uv=False)[-1])
        lam2s.append(info['second_multiplier'])
        assert info['null_residual'] < 1e-7, \
            'the residual should stay clean -- that is the whole point'

    assert lam2s[0] > 0.98 and lam2s[1] > 0.999, \
        'the slow node did not produce a near-unit multiplier (%s), so ' \
        'nothing is being tested' % lam2s
    decay = smins[0] / smins[1]
    assert 30 < decay < 300, \
        'sigma_min of the bordered system fell %.1fx for 100x the time ' \
        'constant (%.3e -> %.3e); the degradation is supposed to track ' \
        'T/tau' % (decay, smins[0], smins[1])


def test_the_ppv_says_when_a_second_multiplier_is_near_one():
    """The detector, because no residual can report this.

    Deflated power iteration against the eigenvectors `ppv()` already has:
    `u` spans the unit mode and `v` is its left partner, so the projection
    removes it exactly and what is left converges to `|lambda_2|`.
    Recovered to six digits against a dense eigendecomposition.

    ⚠ THE WARNING IS ABOUT MORE THAN CONDITIONING. The phase equation
    treats the oscillator's frequency response as INSTANTANEOUS, so a slow
    node that filters a nearby device's noise is invisible to it and phase
    noise comes out OVER-ESTIMATED. Better extraction does not fix that —
    Lai is explicit that the PPV "can be extracted correctly"
    and the analysis is "still inaccurate". So the result is an upper
    bound, and the warning says so rather than implying a tolerance would
    help.
    """
    for tt, expect in ((None, False), (1e0, False), (1e2, True)):
        cir, pss = _solve_slow(tt)
        fp = pss.factored_period()
        n = fp.width
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        true2 = np.sort(np.abs(np.linalg.eigvals(M)))[::-1][1]
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            _v, info = pss.ppv()
        got = info['second_multiplier']
        assert abs(got - true2) < 1e-4 * max(true2, 1e-3), \
            'the deflated power iteration reports |lambda_2| = %.6f ' \
            'against %.6f from a dense eigendecomposition' % (got, true2)
        fired = any('SECOND Floquet multiplier' in str(w.message)
                    for w in rec)
        assert fired is expect, \
            'tau/T=%r: |lambda_2| = %.6f and warned=%s' % (tt, got, fired)


@pytest.mark.parametrize('tau_over_T', [10.0, 100.0])
def test_the_oscillator_Q_is_recovered_from_the_second_multiplier(tau_over_T):
    """⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS, and it is already computed.

    An amplitude perturbation decays to `|lambda_2|` of its size each
    cycle, so the cycles needed to fall below a threshold IS the
    oscillator's Q: `Q = log(threshold)/log|lambda_2|`. The usual
    definitions do not apply to an autonomous circuit — `f_r/df` presumes a
    Bode plot of a BIBO-stable linear system, stored/dissipated presumes
    damping a self-sustaining oscillator does not have, and it is NOT the
    Q of the resonator inside it.

    ⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS MODULE WARNS
    ABOUT. "High Q", "a second multiplier near 1", "slow amplitude
    restoration" and "a long time constant" are four vocabularies for one
    thing — which is why the same circuits defeat the phase row, the
    eigen-split, the probe's continuation and the PPV's
    instantaneous-response assumption.

    ⚠ THE TEST IS A CONSTRUCTION CHECK, WHICH IS WHY IT IS WORTH HAVING.
    The slow node is built with a KNOWN `tau/T`, and its multiplier is
    `exp(-T/tau)`, so `Q` must come back as `tau/T` itself. Measured 10.00
    and 99.99 for 10 and 100 — an independent quantity reproducing an input
    the computation never sees.
    """
    _cir, pss = _solve_slow(tau_over_T)
    with quiet(AccuracyWarning):
        _v, info = pss.ppv()
    Q = info['Q']
    assert abs(Q - tau_over_T) < 0.02 * tau_over_T, \
        'Q came back %.3f for a node built with tau/T = %g; Q is defined ' \
        'as cycles-to-1/e and the multiplier is exp(-T/tau), so the two ' \
        'are the same number' % (Q, tau_over_T)
    assert info['second_multiplier'] > 0.9, \
        'the slow node did not produce a near-unit multiplier'


def test_a_fast_oscillator_has_a_small_Q():
    """The control: plain van der Pol restores amplitude within a cycle.

    `|lambda_2| = 8.5e-04`, so Q ~ 0.14 — the perturbation is gone before
    the cycle ends, which is what "low Q" means for an oscillator and is
    the opposite corner from the case that breaks every method here.
    """
    _cir, pss = _solve_slow(None)
    with quiet():
        _v, info = pss.ppv()
    assert info['Q'] < 1.0, \
        'plain van der Pol reports Q = %.3f; it restores amplitude in ' \
        'well under a cycle' % info['Q']


def _vdp_scaled(cval, lval, period, npts=400):
    """van der Pol with reactances that are NOT unity — see the test below."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=cval)
    c['L'] = L('v', gnd, L=lval)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=period, timestep=period / npts,
                  x0=np.array([2.0, 0.0]), maxiterations=80)
    assert pss.converged
    v, info = pss.ppv()
    irn = pss.irefnode
    x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
    x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
    Cf = np.asarray(c.C(x0f))
    Cm = np.delete(np.delete(Cf, irn, 0), irn, 1).astype(float)
    return c, pss, v, info, Cm


def test_the_ppv_fixture_is_blind_to_a_missing_C_from_one_direction():
    """⚠ THIS TESTS THE TEST, and it found a real hole in the fixture.

    `v·δ` (right) and `vᵀCδ` (the transcription of Demir's Remark 3.1 that
    cost 7% before it was measured out) are distinguished by
    `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has` —
    but only because it perturbs in four RANDOM directions. Van der Pol as
    shipped has `C = diag(1, −1)`, so along the capacitor node the two
    formulations are **numerically identical**: `v·e₀ = vᵀCe₀` exactly,
    ratio 1.0000. Any single-direction probe at the capacitor — which is
    what two earlier gates in this campaign actually did — cannot see the
    difference at all.

    ⚠ "VERIFIED TO 1e-15" SAYS NOTHING ABOUT WHICH ERRORS A FIXTURE CAN
    SEE. A unit reactance makes `C` the identity up to a sign, and an
    implementation that drops the `C` weighting then passes every check
    exactly. The repair is to run the same circuit at a reactance that is
    not 1: with `c = 2` the same blind direction separates the two
    formulations by exactly the capacitance.

    Asserted on both fixtures on purpose — the blindness is recorded as a
    measured property of the shipped one, not as a hypothesis about it, so
    that a future test written against van der Pol knows what it is
    choosing when it perturbs at the capacitor.
    """
    _c1, _p1, v1, _i1, Cm1 = _vdp_scaled(1.0, 1.0, 6.6634)
    m = 2
    e0 = np.array([1.0, 0.0])
    plain1 = float(v1[:m] @ e0)
    weighted1 = float(v1[:m] @ (Cm1 @ e0))
    assert abs(weighted1 / plain1 - 1.0) < 1e-12, \
        'the shipped fixture was expected to be BLIND here (ratio 1); it ' \
        'now reads %.6f, so C is no longer diag(1,-1) and this test has ' \
        'stopped describing the fixture' % (weighted1 / plain1)

    _c2, _p2, v2, _i2, Cm2 = _vdp_scaled(2.0, 3.0, 6.6634 * np.sqrt(6.0))
    plain2 = float(v2[:m] @ e0)
    weighted2 = float(v2[:m] @ (Cm2 @ e0))
    assert abs(weighted2 / plain2 - 2.0) < 1e-9, \
        'ratio %.6f; at c = 2 the missing-C error must show as exactly ' \
        'the capacitance from the SAME direction the unit fixture cannot ' \
        'see it from' % (weighted2 / plain2)
    ## and the normalisation itself must survive the rescaling -- the
    ## property under test is the fixture's discriminating power, not a
    ## claim that a scaled circuit is solved differently
    assert abs(float(v2[:m] @ _i2['xdot']) - 1.0) < 1e-9


def test_the_ppv_gate_probes_the_direction_of_maximum_sensitivity():
    """⚠ THE DESIGNED PROBE, replacing four random directions and hope.

    `v·δ` and `vᵀCδ` agree exactly when `vᵀ(C − I)δ = 0`, so the set of
    directions blind to a dropped `C` is a HYPERPLANE with a known normal:

        blind  ⟺  δ ⊥ (C − I)ᵀv

    which makes `(C − I)ᵀv` itself the direction of maximum sensitivity,
    available for one matrix-vector product from quantities already in
    hand. Van der Pol's `C = diag(1, −1)` gives `(C − I)ᵀv = (0, −2v₁)`, so
    `e₀` — the capacitor node — is orthogonal to it EXACTLY. That is why
    the companion test measures a ratio of 1.0000 there: not a near miss,
    an exact one.

    ⚠ ALONG THE DESIGNED PROBE THE TWO HYPOTHESES PREDICT OPPOSITE SIGNS,
    so no tolerance is needed to separate them — the circuit picks one:

        npts   measured      v·δ (right)   vᵀCδ (wrong)   ratio
         400   −5.4215e-01   −5.4131e-01   +5.4131e-01    0.998443
         800   −5.4153e-01   −5.4110e-01   +5.4110e-01    0.999213

    The magnitude converges at O(h) as a bonus; the SIGN alone already
    decides it. Four random directions scored 0.43/0.78/0.81/0.96 on this
    fixture — they work, and the spread is exactly the luck this removes.

    ⚠ THIS IS A TARGETED PROBE AND ITS OPTIMALITY IS ABOUT ONE ERROR. It
    maximises sensitivity to a dropped `C` weighting and says nothing
    about any other defect; the random-direction gate stays because it is
    not aimed at a hypothesis. Recording which errors a check can see is
    the whole point of §D shape 0c, and that applies to this one too.

    ⚠ AND THE STRUCTURAL RULE IS WORTH MORE THAN THIS CIRCUIT. `C` is ZERO
    on algebraic rows, so `C − I = −I` there and the discrepancy is
    maximal: in an MNA-shaped system the rows a DAE solver already treats
    specially are the BEST probes, and the capacitor nodes everyone
    reaches for first are the blind ones.
    """
    from pycircuit.circuit.transient import Transient
    out = []
    for npts in (400, 800):
        cir, pss, v, info = _vdp_ppv(npts)
        m = cir.n - 1
        irn = pss.irefnode
        T = pss.period
        xdot = info['xdot']
        x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        Cm = np.delete(np.delete(np.asarray(cir.C(x0f)), irn, 0),
                       irn, 1).astype(float)

        normal = (Cm - np.eye(m)).T @ v[:m]
        assert abs(float(normal @ np.eye(m)[0])) < 1e-12 * max(
            np.linalg.norm(normal), 1e-300), \
            'the capacitor axis is no longer exactly in the blind ' \
            'hyperplane, so this fixture has changed shape'
        d = normal / np.linalg.norm(normal)

        def integrate(xi, nper=3, ppp=2000):
            tran = Transient(cir, toolkit=circuit.numeric, reltol=1e-9,
                             iabstol=1e-13, vabstol=1e-11)
            with quiet():
                res = tran.solve(refnode=gnd, tend=nper * T,
                                 timestep=T / ppp, x0=xi)
            return np.asarray(res.x, dtype=float)[:, -1]

        eps = 1e-5
        ref = integrate(x0f)
        dr = np.concatenate((d[:irn], np.zeros(1), d[irn:]))
        dx = np.delete(integrate(x0f + eps * dr) - ref, irn)
        meas = float(dx @ xdot) / float(xdot @ xdot) / eps
        p_ok = float(v[:m] @ d)
        p_bug = float(v[:m] @ (Cm @ d))

        assert np.sign(p_ok) != np.sign(p_bug), \
            'the designed probe no longer separates the two formulations ' \
            'by sign; it is not the maximum-sensitivity direction'
        assert np.sign(meas) == np.sign(p_ok), \
            'the CIRCUIT chose the vTCd formulation (measured %+.4e, ' \
            'v.d %+.4e); the PPV normalisation is wrong' % (meas, p_ok)
        out.append(abs(p_ok / meas - 1.0))

    assert out[0] < 5e-3, 'coarse grid off by %.3e' % out[0]
    assert out[1] < 0.75 * out[0], \
        'the designed probe does not converge (%.3e -> %.3e); a residue ' \
        'that does not shrink is a defect, not discretisation'\
        % (out[0], out[1])


def _vdp_at_Q(Q, npts=480, mu=None):
    """A van der Pol tuned to a target `Q` — `μ = 1/(2πQ)`.

    ⚠ THE RECIPE IS MEASURED, NOT ASSUMED. `λ₂ ≈ exp(−μT)` with `T ≈ 2π`
    gives `Q = 1/(μT) = 1/(2πμ)`, and against this solver: predicted
    3.183/7.958/15.92/63.66 against measured 3.182/7.959/15.92/63.67 at
    `μ` = 0.05/0.02/0.01/0.0025 — four digits.

    ⚠ THE PERIOD SEED MATTERS AT SMALL `μ`. `2π/√(1−μ²/4)` and
    `reltol = 1e-12`, or the shooting solve becomes the thing under test
    rather than the instrument measuring it.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q) if mu is None else mu
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                  maxiterations=150)
    assert pss.converged, 'mu = %r did not converge' % mu
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_reported_Q_amplifies_its_own_lambda2_error_by_Q():
    """⚠⚠ `info['Q']` CARRIES A `Q`-FOLD AMPLIFIED ERROR, AND NOTHING SAID SO.

    `Q = −1/ln λ₂` differentiates to

        (dQ/Q) / (dλ₂/λ₂)  =  −1/ln λ₂  =  Q

    so **the relative error in `Q` is `Q` times the relative error in
    `λ₂`**. That makes the identity behind §0 of the roadmap do double
    duty: it is what ties the seven `λ₂` failure modes together, AND it
    multiplies every `λ₂` error by `Q` on the way out.

    ⚠ MEASURED END TO END IN THIS SOLVER, not just in the algebra —
    `λ₂` against the finest grid, and the ratio of the two relative
    errors:

        Q_target   npts   rel err λ₂   rel err Q    ratio
          3.18      120    1.04e-03     3.31e-03      3.2
         15.92      120    5.75e-04     9.23e-03     16.1
         63.66      120    4.89e-04     3.21e-02     65.7
         63.66      240    6.29e-05     4.02e-03     63.9
         63.66      480    7.56e-06     4.81e-04     63.7

    ⚠ SO THE EXPOSURE IS A RESOLUTION REQUIREMENT THAT SCALES WITH `Q`,
    not a fixed accuracy. 120 points/period gives `Q` to 0.3% at `Q = 3`
    and to only 3.2% at `Q = 64`; to report `Q` to 1% at `Q = 100` needs
    `λ₂` to 1e-4 relative, and at `Q = 1000` to 1e-5.

    ⚠ AND THE GOOD NEWS IS THAT IT IS A RESOLUTION PROBLEM RATHER THAN A
    BIAS. Gear-2's `λ₂` converges here at better than second order
    (4.89e-04 → 6.29e-05 → 7.56e-06, ~8× per doubling), so the requirement
    is payable. A method that biased `λ₂` at fixed order — backward Euler
    does — would not have that escape, and the amplification would turn a
    5.6e-2 bias into 85% at `Q = 100`.

    ⚠ THIS TEST DOES NOT ASSERT ACCURACY. It asserts the SENSITIVITY, so
    that the cost of reading `Q` is pinned rather than discovered. A
    future change that broke the identity would show up here as a ratio
    that is no longer `Q`.
    """
    for Q in (3.183, 15.915, 63.662):
        rows = []
        for npts in (120, 480):
            _cir, pss, _pac = _vdp_at_Q(Q, npts=npts)
            _v, info = pss.ppv()
            rows.append((info['second_multiplier'], info['Q']))
        (l_c, q_c), (l_f, q_f) = rows
        el = abs(l_c / l_f - 1.0)
        eq = abs(q_c / q_f - 1.0)
        assert el > 0, 'lambda2 identical at two grids; nothing to amplify'
        ratio = eq / el
        assert abs(ratio / q_f - 1.0) < 0.10, \
            'at Q = %.3f the error amplification is %.2f, not Q. Either ' \
            'info["Q"] is no longer -1/log(lambda2) or the identity ' \
            'Q = log(threshold)/log|lambda2| has been broken' % (q_f, ratio)
        ## and the recipe itself, so the fixture cannot drift
        assert abs(q_f / Q - 1.0) < 0.02, \
            'mu = 1/(2 pi Q) gave Q = %.4f against %.4f requested; the ' \
            'fixture recipe no longer holds on this solver' % (q_f, Q)


def _vdp_with_parasitic(Q, tau_over_T, npts=480):
    """A high-`Q` van der Pol with one parasitic RC node at a chosen `tau/T`.

    Lets a parasitic multiplier be swept THROUGH the oscillatory one:
    `lambda_p = exp(-T/tau_p)` equals `lambda_2 = exp(-1/Q)` exactly when
    `tau_p/T = Q`. A resonance between a designed quantity and an
    incidental one.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir.add_node('w')
    rbig = 1e6
    cir['Rs'] = R('v', 'w', r=rbig)
    cir['Cs'] = C('w', gnd, c=tau_over_T * T / rbig)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=150)
    assert pss.converged, 'Q=%r tau/T=%r' % (Q, tau_over_T)
    return cir, pss


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_a_resonant_parasitic_breaks_eigenvectors_and_not_the_ppv():
    """⚠⚠ THE BORDERED SOLVE EARNS ITS COST HERE, MEASURED.

    A parasitic multiplier `λ_p = exp(−T/τ_p)` equals `λ₂ = exp(−1/Q)`
    exactly when `τ_p/T = Q` — **a resonance between a quantity the
    designer chooses and one they do not**. At `Q = 16`, `τ_p/T = 16`
    gives `λ₂ = 0.939432` against `λ_p = 0.939410`: degenerate to five
    digits.

    Swept through it:

        τ_p/T    PPV drift    λ₂ eigenvector drift   cond(V)
          1      2.1e-08      0                       4.7
          8      2.6e-08      0.0251                  4.9
         15      3.0e-08      0.0254                  7.5
         16      3.3e-08      0.0254                 25.7   ← resonance
         17      3.0e-08      0.99967                 4.5   ← 90° SWAP
         32      3.0e-08      0.99967                10.8

    ⚠ **The λ₂ eigenvector swaps by 90° across the crossing and the PPV
    does not move at all** — flat at 3e-08 through resonance, with the
    border and null residuals unchanged at 1.8e-13 / 7.0e-13.

    ⚠⚠ AND THE REASON IS THE ONE THAT MAKES THE DISTINCTION USEFUL: the
    degeneracy is between `λ₂` and `λ_p`, **neither of which is 1**. The
    PPV is the left null vector of `I − M`, i.e. the `λ₁ = 1` object,
    which stays SIMPLE throughout. So a degeneracy among the *other*
    multipliers cannot touch it.

    ⚠ THAT IS NOT DEMIR & ROYCHOWDHURY'S OBJECTION, AND CONFLATING THEM
    WOULD MIS-SCOPE BOTH. Theirs is `λ₂ → λ₁ = 1` — the HIGH-Q case,
    where the phase mode itself becomes indistinguishable, which is what
    `PPV_SECOND_MULTIPLIER_WARN` guards. This is a different degeneracy
    with a different victim: it breaks any eigen-based extraction of `λ₂`
    and leaves the phase mode alone.

    Two conclusions, both worth having separately: an eigendecomposition
    route to the PPV would be unreliable here for a reason the *residual*
    cannot see, and the bordered solve is immune by construction rather
    than by luck.
    """
    ref_v = None
    ref_u2 = None
    out = []
    for tau in (1.0, 16.0, 32.0):
        _cir, pss = _vdp_with_parasitic(16.0, tau)
        fp = pss.factored_period()
        n = fp.width
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        w, V = np.linalg.eig(M)
        order = np.argsort(-np.abs(w))
        w, V = w[order], V[:, order]
        v, info = pss.ppv()
        vv = np.real(v) / np.linalg.norm(v)
        u2 = np.real(V[:, 1])
        u2 = u2 / np.linalg.norm(u2)
        if ref_v is None:
            ref_v, ref_u2 = vv.copy(), u2.copy()

        def sin_angle(a, b):
            c = abs(float(a @ b))
            return np.sqrt(max(1.0 - c ** 2, 0.0))

        out.append((tau, sin_angle(vv, ref_v), sin_angle(u2, ref_u2),
                    abs(w[1] - w[2]), info['null_residual'],
                    abs(info['border_residual'])))

    for tau, dv, du, gap, nres, bres in out:
        assert dv < 1e-6, \
            'tau/T = %g moved the PPV by sin = %.3e; the phase mode is ' \
            'supposed to be untouched by a lambda_2/lambda_p degeneracy ' \
            'because neither of them is 1' % (tau, dv)
        assert nres < 1e-9 and bres < 1e-9, \
            'tau/T = %g degraded the bordered solve (null %.2e, border ' \
            '%.2e); it is supposed to be immune by construction' \
            % (tau, nres, bres)
    ## the resonance really is a near-degeneracy, and it really does move
    ## the eigenvector -- otherwise the test above proves nothing
    gaps = {t: g for t, _dv, _du, g, _n, _b in out}
    assert gaps[16.0] < 1e-3, \
        'tau/T = 16 no longer puts lambda_p on top of lambda_2 (gap ' \
        '%.3e), so this fixture has stopped testing the resonance' \
        % gaps[16.0]
    assert gaps[1.0] > 0.1 and gaps[32.0] > 0.01
    assert out[2][2] > 0.5, \
        'the lambda_2 eigenvector no longer swaps across the crossing ' \
        '(sin = %.3e); if eigen extraction has become stable here the ' \
        'contrast this test rests on is gone' % out[2][2]


def _high_q_with_bulk(Q=60.0, nbulk=10, lam_lo=0.05, lam_hi=0.35,
                      rbig=1e6, psd=1e-6, npts=480):
    """A high-`Q` oscillator with a genuine DAMPED BULK — `m = 12`.

    ⚠ THE FIXTURE THIS FILE SPENT A LONG TIME NOT HAVING. Every gate here
    was written against van der Pol at `μ = 1`: `λ₂ = 8.6e-4` and `m = 2`.
    Three separate questions ended in "we need a different fixture, not a
    better probe" — the Ritz gate (no bulk), the PPV gate's `1/√m`
    protection, and the transverse-kick degeneracy.

    Design constraints, each of which is a lesson:

    * `μ = 1/(2πQ)` sets the core — verified to four digits;
    * the bulk must be **FAST**. Slow RC nodes add *more* near-unit
      multipliers, which is the opposite of a bulk and would break the
      one-dimensional null space `ppv` and `oscillator_covariance` rest
      on. `λ_p` is placed log-uniformly in `[0.05, 0.35]`;
    * `τ_p/T ≪ Q` avoids the resonance `τ_p/T = Q` where `λ_p` collides
      with `λ₂`;
    * coupling through `rbig` so the branches do not load the orbit.

    MEASURED: `m = 12`, `n = 24`, `Q = 60.24`, `λ₂ = 0.9835`, bulk
    0.0500–0.3500, `cond(V) = 92`, converges in ~4 s.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    for i, lp in enumerate(np.exp(np.linspace(np.log(lam_lo),
                                              np.log(lam_hi), nbulk))):
        nm = 'p%d' % i
        cir.add_node(nm)
        cir['Rp%d' % i] = R('v', nm, r=rbig)
        cir['Cp%d' % i] = C(nm, gnd, c=(-T / np.log(lp)) / rbig)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=200)
    assert pss.converged, 'the high-Q bulk fixture did not converge'
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_high_q_bulk_fixture_is_what_it_claims():
    """The fixture's own regression test — it is useless if it drifts.

    Four properties, each of which a plausible edit would break: the core
    sets `Q`, the bulk is FAST and spread, `λ₂` is the oscillator's mode
    and not a parasitic, and the eigenvector conditioning stays modest.
    """
    _cir, pss, _pac = _high_q_with_bulk()
    m = pss.cir.n - 1
    assert m == 12, 'm = %d; the 1/sqrt(m) exposure needs a dozen states' % m
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    ev = np.sort_complex(np.linalg.eigvals(M))[::-1]
    _v, info = pss.ppv()
    assert abs(info['Q'] / 60.0 - 1.0) < 0.02, \
        'Q = %.3f, not 60; mu = 1/(2 pi Q) no longer sets the core' \
        % info['Q']
    ## the bulk: ten modes, fast, spread, and well clear of lambda_2
    bulk = np.abs(ev[2:12])
    assert bulk.max() < 0.4, \
        'the bulk reaches %.4f; slow parasitics are MORE near-unit modes, ' \
        'not a bulk, and would break the 1-D null space' % bulk.max()
    assert bulk.min() > 0.02 and bulk.max() / bulk.min() > 3.0, \
        'the bulk is not spread (%.4f..%.4f)' % (bulk.min(), bulk.max())
    assert abs(ev[1]) > 0.95, \
        'lambda_2 = %.4f is no longer the oscillator amplitude mode' \
        % abs(ev[1])
    ## and the conditioning the Ritz filter depends on
    cond = float(np.linalg.cond(np.linalg.eig(M)[1]))
    assert cond < 1e3, \
        'cond(V) = %.3e; above ~1e4 the |lam-1| filter in ppv() stops ' \
        'resolving the phase mode and would select it as lambda_2' % cond


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_ppv_physical_gate_cannot_verify_the_ppv_at_high_q():
    """⚠⚠ THE ONE GATE THAT DOES BREAK, AND NO CHEAP REPAIR WORKS.

    On the `m = 12`, `Q = 60` fixture the physical gate reads **24% wrong
    and does not converge**:

        npts   random dirs   |b/a| = 1 constructed
         240   2.66e-01      8.46e-02
         480   2.42e-01      7.71e-02
        ratio  1.097         1.097          (2.0 would be O(h))

    ⚠ IT IS NOT DISCRETISATION. A ratio of 1.097 across a grid doubling
    says `npts` cannot touch it. It is transverse contamination, and at
    `λ₂ = 0.9835` essentially none of it decays in one period — the honest
    cost of waiting it out is `4.6·Q ≈ 277` periods.

    ⚠ THE ERROR TRACKS THE DIRECTION, NOT THE WAIT. Per-direction at
    `npts = 480`: `|b/a|` = 31.3 → 98% error, 5.0 → 29%, 1.13 → 4.0%. And
    `1/√m` makes a random draw in 12 dimensions strongly transverse, which
    is exactly the protection van der Pol had and this does not.
    Constructing `|b/a| = 1` cuts 24% to 8% — real, and not a fix.

    ⚠⚠ AND EXTRAPOLATION FAILS HERE TOO, FOR A NEW REASON. Median over
    four directions: raw n=1 **16.4%**, raw n=3 46.6%, Aitken **85.9%**,
    `λ₂`-extrapolation **110%**. At `m = 2` extrapolation failed because
    the residual was O(h) rather than the geometric mode; here it fails
    because there are **eleven** transverse modes and Aitken models
    exactly one. **The bulk that makes the fixture realistic is what
    defeats the repair.**

    ⚠ SO WHAT THIS ESTABLISHES IS A LIMIT ON THE GATE, NOT A DEFECT IN THE
    PPV. The bordered solve's residuals on this fixture are 2.0e-13 and
    5.1e-14, unchanged from `m = 2`. The PPV may well be right; **this
    experiment cannot say so at high Q**, and no cheap variant of it can.
    Recorded so that "the PPV is gated physically" is not read as covering
    the regime §0 says matters.
    """
    _cir, pss, _pac = _high_q_with_bulk()
    _v, info = pss.ppv()
    ## the bordered solve itself is untroubled -- that is the point
    assert info['null_residual'] < 1e-9, \
        'null residual %.2e; the SOLVE is supposed to be fine here' \
        % info['null_residual']
    assert abs(info['border_residual']) < 1e-9
    ## and a random direction really is strongly transverse at m = 12
    m = pss.cir.n - 1
    xh = np.asarray(info['xdot'], dtype=float)
    xh = xh / np.linalg.norm(xh)
    rng = np.random.default_rng(0)
    ratios = []
    for _ in range(4):
        d = rng.standard_normal(m)
        d /= np.linalg.norm(d)
        tan = abs(float(d @ xh))
        ratios.append(np.sqrt(max(1.0 - tan ** 2, 0.0)) / max(tan, 1e-300))
    assert min(ratios) > 1.0, \
        'a random direction at m = %d should be mostly transverse ' \
        '(1/sqrt(m) = %.2f tangential); got |b/a| min %.2f' \
        % (m, 1.0 / np.sqrt(m), min(ratios))


def _asym_lossy_osc(Q, idc=0.0, npts=480, seedT=None,
                    arel=0.25, rrel=0.2):
    """An oscillator with `∫v₀ dt ≠ 0`, at any `Q`.

    ⚠ BOTH PERTURBATIONS SCALE WITH `μ`, and that is not cosmetic. At
    `Q = 30`, `μ = 5.3e-3`; a fixed `rs = 0.02` contributes an effective
    conductance **four times** the negative resistance and the oscillator
    does not start, while a fixed `a = 0.25` swamps the van der Pol term
    entirely. Written as `a = 0.25μ`, `rs = 0.2μ` this reduces at `μ = 1`
    to exactly the A4d fixture, and stays a *perturbed* oscillator as `Q`
    rises.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0)
                       + arel * mu * (u ** 2 - 2.0))
    cir.add_node('x')
    cir['L'] = L('v', 'x', L=1.0)
    cir['Rs'] = R('x', gnd, r=rrel * mu)
    if idc:
        cir['I'] = IS('v', gnd, i=idc)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=seedT or 2.0 * np.pi, timestep=(seedT or 2.0
                  * np.pi) / npts, x0=x0, maxiterations=250)
    assert pss.converged, 'Q=%r idc=%r did not converge' % (Q, idc)
    return cir, pss


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_ppv_predicts_a_frequency_shift_at_high_q():
    """⚠⚠ THE PPV GATE THAT WORKS AT HIGH Q — because it has no transient.

    The state-kick gate cannot verify the PPV at `Q = 60`: its observable
    is contaminated by transverse modes that need `4.6·Q ≈ 277` periods to
    decay. This one perturbs a SOURCE and measures the resulting PERIOD:

        ΔT / δi  =  ∫₀ᵀ v₀(t) dt

    ⚠ THE NEW PERIOD IS A PROPERTY OF THE CONVERGED ORBIT, so there is no
    transient for a slow mode to contaminate. The measurement re-solves
    the PSS and never touches `v`, the monodromy or the bordered solve.

    MEASURED at `Q = 75`:

        idc     ΔT/ulp     ratio to prediction
        1e-06        3.2   0.967611
        1e-04      328.5   0.998574
        1e-02    32896.7   1.000006      ← six parts per million
        1e-01   328992.4   1.000084
        3e-01   987106.3   1.000215

    ⚠ A TWO-SIDED OPTIMUM, AND BOTH SIDES ARE UNDERSTOOD. Below ~300 ulp
    the period difference is round-off — at `idc = 1e-6` it is 3.2 ulp,
    and there the error grows with `npts` because more steps accumulate
    more of it, which is the tell. Above `1e-2` the second-order term in
    `δi` takes over, linear in `idc` as it should be (8.4e-5 → 2.15e-4 for
    a 3× injection). A single-amplitude assertion would have landed on one
    side or the other and read as a defect.

    ⚠⚠ AND THIS GATE IS THE SAME FUNCTIONAL AS A4d's `Γ`, which is why it
    exists at all here: a DC current IS a zero-frequency coloured source,
    so both contract `∫v dt`. It therefore validates
    `colour_projection`'s kernel — which had no external gate — and it
    inherits A4d's preconditions exactly. On a lossless symmetric tank
    both the prediction and the measurement are zero (measured `4e-11` and
    `1e-15`), which is a confirmation of the structural identity and
    useless as a gate.
    """
    Q = 60.0
    _c0, p0 = _asym_lossy_osc(Q)
    T0 = float(p0.period)
    m = p0.cir.n - 1
    ## ⚠ the raw pair, whose integral carries the same-grid identity this
    ## gate pins; the consistent `samples` has a ~1e-5 |v| floor on its
    ## mean and the true mean here is 4.6e-9 (see `_raw_pair_integrals`)
    ints, info = _raw_pair_integrals(p0, _c0)
    pred = ints[0]
    assert info['Q'] > 40.0, \
        'Q = %.2f; this gate exists to run in the regime the state-kick ' \
        'gate cannot reach' % info['Q']
    assert abs(pred) > 1e-9, \
        'int v0 dt = %.3e is at the structural zero; an asymmetric LOSSY ' \
        'tank is required or there is nothing to measure' % pred

    ratios = []
    for idc in (1e-4, 1e-2, 1e-1):
        _c1, p1 = _asym_lossy_osc(Q, idc=idc, seedT=T0)
        ratios.append(((float(p1.period) - T0) / idc) / pred)
    ## the sweet spot: off the round-off floor, below the quadratic term
    assert abs(ratios[1] - 1.0) < 1e-4, \
        'at idc = 1e-2 the PPV predicts dT/di to %.3e relative; this is ' \
        'the one measurement that reaches high Q with no transient in it' \
        % abs(ratios[1] - 1.0)
    ## and both edges behave as their mechanisms say
    assert abs(ratios[0] - 1.0) > abs(ratios[1] - 1.0), \
        'the small injection is no longer round-off limited, so the ' \
        'two-sided structure this test documents has changed'
    assert abs(ratios[2] - 1.0) > abs(ratios[1] - 1.0), \
        'the large injection no longer shows the second-order term'


def test_the_frequency_shift_gate_converges_at_first_order():
    """`O(h)` on the low-`Q` fixture, where round-off is far away.

    At `μ = 1` the signal is `∫v₀ dt ≈ 9.9e-4` — five orders above the ulp
    floor — so the residual is pure discretisation and must halve with the
    grid. Measured 9.92e-05 → 4.58e-05, **ratio 2.17**.

    Asserted separately from the high-`Q` test because they establish
    different things: this one that the gate is a *converging* measurement
    of the right quantity, that one that it still works where the
    alternative does not.
    """
    errs = []
    for npts in (480, 960):
        _c0, p0 = _asym_lossy_osc(1.0, npts=npts)
        T0 = float(p0.period)
        ## the raw pair -- see `_raw_pair_integrals`
        ints, _info = _raw_pair_integrals(p0, _c0)
        pred = ints[0]
        _c1, p1 = _asym_lossy_osc(1.0, idc=1e-6, npts=npts, seedT=T0)
        errs.append(abs(((float(p1.period) - T0) / 1e-6) / pred - 1.0))
    assert errs[0] < 3e-4, 'coarse grid off by %.3e' % errs[0]
    assert errs[1] < 0.6 * errs[0], \
        'the gate does not converge (%.3e -> %.3e); a residual that does ' \
        'not shrink is a defect, not discretisation' % (errs[0], errs[1])


def test_the_ppv_is_invariant_to_the_newtons_inner_solver():
    """⚠ THE PROPERTY THAT MAKES B6 UNNECESSARY, pinned so it cannot regress.

    García, Romero & Acha get the Floquet multipliers from the Hessenberg
    matrix their shooting-Newton's GMRES already built. We adopted the
    **Ritz values** (`fef3d60`) and kept a small dedicated Arnoldi rather
    than reusing the Newton's basis — and the reason is measured, not
    assumed:

    ⚠ THE DEFAULT NEWTON HAS NO GMRES AT ALL. `solve(matrix_free=False)`
    is the default and factors the Jacobian directly, so there is no
    Hessenberg matrix to reuse. Instrumenting the bulk fixture's solve
    counted **zero** inner GMRES matvecs.

    ⚠ AND WHERE THERE IS ONE, THE SAVING IS 12 MATVECS — measured at 27%
    of `ppv()` but only **~1.5%** of a PSS-plus-`ppv` workflow (Arnoldi
    0.077 s against PSS 4.17 s + `ppv` 0.29 s on the `m = 12` fixture),
    shrinking further on the large circuits where `matrix_free` is
    actually worth using, because the PSS dominates more there. Realising
    it means replacing scipy's `gmres` in the core Newton — the
    highest-risk change available — for that.

    ⚠⚠ AND NOTHING DEPENDS ON THE CHOICE, WHICH IS WHAT THIS TEST HOLDS.
    `factored_period()` re-traverses at the CONVERGED solution, so the
    multipliers, `Q` and the PPV come out **bit-identical** either way.
    That is a real design property and a fragile one: the class docstring
    records that `Jtvec`/`Cvec` are written by neither factored traversal,
    so an analysis reading them after a matrix-free solve would rebuild an
    operator for a different trajectory, silently. The re-traversal is
    what keeps that from mattering here.
    """
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    ## ⚠ THE STAGE PREDICTOR IS PINNED OFF, and that is a statement about what
    ## "bit-identical" can mean here.  The re-traversal IS a pure function of
    ## `x_in` either way -- `_begin_period` resets the predictor's node history
    ## for exactly that reason -- but the two inner solvers converge to `x_in`
    ## values that differ in the last ulp, and a predictor's extrapolation
    ## carries that into the traversal instead of contracting it: measured,
    ## lambda_2 agrees to 7e-11 relative rather than to the bit.  The design
    ## property this test names is about the SOLVER, so it is measured on a
    ## fixed seed; the tolerance-level claim with the predictor on is asserted
    ## at the end.
    prev_pred = Transient.stage_predictor
    Transient.stage_predictor = 'off'
    out = {}
    for mf in (False, True):
        cir = _vdp_with_noise(1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        with quiet(AccuracyWarning):
            ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
            pss.solve(period=6.6634, timestep=6.6634 / 239, x0=x0,
                      maxiterations=60, matrix_free=mf)
        assert pss.converged
        with quiet():
            v, info = pss.ppv()
        out[mf] = (np.asarray(v).copy(), info['second_multiplier'],
                   info['Q'], info['null_residual'])

    v0, l0, q0, r0 = out[False]
    v1, l1, q1, r1 = out[True]
    assert l0 == l1, \
        'lambda_2 differs between the dense and matrix-free Newton ' \
        '(%.12g vs %.12g); factored_period() is supposed to re-traverse ' \
        'at the converged solution, so the inner solver cannot matter' \
        % (l0, l1)
    assert q0 == q1
    assert np.array_equal(v0, v1), \
        'the PPV differs between inner solvers; the largest component ' \
        'gap is %.3e' % float(np.max(np.abs(v0 - v1)))
    assert r0 < 1e-9 and r1 < 1e-9

    ## and with the predictor on, where the last ulp of `x_in` is carried
    ## rather than contracted: still the same answer, to 1e-9 relative
    Transient.stage_predictor = 'on'
    try:
        pred = {}
        for mf in (False, True):
            cir = _vdp_with_noise(1e-6)
            pss = PSS(cir, method='gear', reltol=1e-12)
            x0 = np.zeros(cir.n - 1)
            x0[0] = 2.0
            with quiet(AccuracyWarning):
                pss.solve(period=6.6634, timestep=6.6634 / 239, x0=x0,
                          maxiterations=60, matrix_free=mf)
                v, info = pss.ppv()
            pred[mf] = (np.asarray(v).copy(), info['second_multiplier'])
    finally:
        Transient.stage_predictor = prev_pred
    assert abs(pred[True][1] - pred[False][1]) < 1e-9 * abs(pred[False][1]), \
        (pred[True][1], pred[False][1])
    assert np.max(np.abs(pred[True][0] - pred[False][0])) < \
        1e-9 * max(float(np.max(np.abs(pred[False][0]))), 1e-30)
    ## and the predictor did not move the answer either
    assert abs(pred[False][1] - l0) < 1e-9 * abs(l0), (pred[False][1], l0)


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_ppv_is_unchanged_by_the_gmres_swap():
    """The bordered solves moved off scipy; the answers must not move.

    `ppv()`'s two augmented solves now use `_arnoldi_gmres` and judge
    convergence by relative residual rather than by a status code. That is
    a change to *how failure is decided*, not to the mathematics, so the
    values are pinned against what the scipy path produced.
    """
    ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
    _cir, pss, _pac = _vdp_at_Q(16.0, npts=479)
    _v, info = pss.ppv()
    assert abs(info['second_multiplier'] - 0.9394257319) < 1e-9, \
        'lambda_2 = %.10f against the 0.9394257319 the scipy path gave' \
        % info['second_multiplier']
    assert abs(info['Q'] - 16.003453) < 1e-4
    assert info['null_residual'] < 1e-11
    assert abs(info['border_residual']) < 1e-11


def _divider_osc(Q=8.0, npts=480, a=0.25, ratio=10.0, idc_node=None, idc=0.0,
                 period=None):
    """The same tank with its loss split into a DIVIDER -- two algebraic nodes.

    ⚠ THIS FIXTURE EXISTS BECAUSE THE SINGLE-RESISTOR ONE CANNOT SETTLE A
    SIGN.  With one series `R`, `|integral v_0|` and `|r integral v_branch|`
    agree to 1.5e-4, so a prediction that matches in magnitude matches for
    EITHER sign and the agreement is no evidence.  Splitting the loss gives
    two algebraic nodes whose fold coefficients into the branch row are
    `(r1 + r2)` and `r2`, so the topology fixes the RATIO of their PPV
    entries -- and `ratio` chooses it.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2 * np.pi * Q)
    rtot = 0.2 * mu
    r2 = rtot / ratio
    r1 = rtot - r2
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0)
                       + a * mu * (u ** 2 - 2.0))
    cir.add_node('x')
    cir.add_node('y')
    cir['L'] = L('v', 'x', L=1.0)
    cir['R1'] = R('x', 'y', r=r1)
    cir['R2'] = R('y', gnd, r=r2)
    if idc_node is not None:
        cir['Idc'] = IS(idc_node, gnd, i=idc)
    T = 2 * np.pi if period is None else float(period)
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=250)
    assert pss.converged, 'divider did not converge'
    return cir, pss, r1, r2


def _reduced_index(cir, pss, name):
    names = [str(nd) for nd in cir.nodes]
    i = names.index(name)
    irn = pss.irefnode
    assert i != irn
    return i if i < irn else i - 1


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_ppv_carries_the_slaved_sensitivity_on_an_algebraic_row():
    """⚠ An algebraic row's PPV entry is SLAVED to the differential ones,
    and it used to be left at zero.

    The PPV entry for a row IS the phase sensitivity to a perturbation
    entering that row, so this is checkable directly: inject a DC current
    at each node and measure `dT/di` against `integral v_j dt`.

    ⚠ THE DIVIDER IS THE FIXTURE AND THAT IS THE POINT.  A single series
    resistor makes `|integral v_0|` and `|r integral v_branch|` agree to
    1.5e-4, so a sign error is invisible there. Splitting the loss so the
    two algebraic nodes fold in with `(r1 + r2)` and `r2` puts their
    entries in a ratio the topology fixes, and the sign has nowhere to
    hide.
    """
    cir, pss, r1, r2 = _divider_osc(ratio=10.0)
    T0 = float(pss.period)
    m = cir.n - 1
    ## ⚠ `samples_eq`: a DC current injection is an EQUATION-ROW input, so
    ## it contracts with `v_1`, not with the `C^T v_1` in `samples`.  And
    ## the RAW pair for the identity this gate pins; the consistent object
    ## is held within its measured absolute floor below.
    ints, info = _raw_pair_integrals(pss, cir)
    S = np.asarray(info['samples_eq'])[:, :m]
    h = np.diff(np.asarray(info['times'], dtype=float))
    n = len(h)
    ints_c = [float((S[:n, j] * h).sum()) for j in range(m)]
    ix = _reduced_index(cir, pss, 'x')
    iy = _reduced_index(cir, pss, 'y')

    ## the STRUCTURAL half: the ratio the topology fixes
    assert abs(ints[ix] / ints[iy] - (r1 + r2) / r2) < 1e-4, \
        'the two algebraic entries must stand in the ratio (r1+r2)/r2 = ' \
        '%.6f; got %.6f' % ((r1 + r2) / r2, ints[ix] / ints[iy])

    ## the MEASURED half: every row, algebraic or not, on ONE convention
    for name in ('v', 'x', 'y'):
        j = _reduced_index(cir, pss, name)
        _c2, p2, _r1, _r2 = _divider_osc(ratio=10.0, idc_node=name,
                                         idc=1e-4, period=T0)
        meas = (float(p2.period) - T0) / 1e-4
        assert abs(meas / ints[j] - 1.0) < 1e-3, \
            'node %s: the PPV predicts dT/di = %+.9e and the measurement ' \
            'gives %+.9e (ratio %+.7f) -- a NEGATIVE ratio is the sign ' \
            'error this fixture exists to catch' % (name, ints[j], meas,
                                                    meas / ints[j])
        ## the consistent object: the same truth, within its floor.  Its
        ## integral on node v read -7.5e-07 / +1.27e-06 / +1.77e-06 at
        ## 480/960/1920 points against an extrapolated +1.9375e-06 --
        ## second order, but the true mean here is 4e-6 |v| and the floor
        ## is ~1e-5 |v|, so at 480 points the SIGN is not resolved.
        assert abs(ints_c[j] - meas) < 5e-6, \
            'node %s: the consistent object is %.2e from the measured ' \
            'dT/di, outside its measured floor' % (name, ints_c[j] - meas)


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_algebraic_fill_is_identified_not_merely_validated():
    """⚠⚠ IDENTIFICATION, WHICH THREE AGREEING REFERENCES ARE NOT.

    The algebraic fill was gated by OUTCOME -- an equivalent circuit, the
    Lyapunov route, and a DC probe all agreed. That leaves open what the
    filled entries ARE. Demir 2000's adjoint, eq (24), settles it:

        C^T(t) dy/dt - G^T(t) y = 0

    ⚠ THE DERIVATIVE IS ON `y` ALONE, NOT ON THE PRODUCT `C^T y` -- the
    contrast Demir flags against his eq (19), and the asymmetry that made a
    from-scratch derivation of this fill come out sign-inverted.

    Row `i` of that system is `(column i of C)^T ydot = (column i of G)^T y`.
    For an ALGEBRAIC state the column of `C` is zero, the left side vanishes,
    and the row degenerates to a POINTWISE CONSTRAINT with no time
    derivative in it at all:

        (column i of G)^T v_1(t) = 0        for every algebraic i, every t

    ⚠⚠ AND THE TEST IS ON `v_1`, WHILE `ppv()` RETURNS `C^T v_1`. On this
    fixture `C = diag(1, 0, -L)`: the INDUCTOR BRANCH ROW CARRIES `-L`, so
    reading the returned vector as `v_1` flips that row. That is exactly why
    the fill's sign had to be flipped against a measurement -- the
    derivation and the measurement were describing different vectors, and
    both were right.

    Measured: the constraint holds at 0.0 EXACTLY once the differential rows
    are divided by their `C` entries, and is violated at O(1) both by the
    raw vector and by the opposite sign. So this discriminates, and it goes
    through none of the three outcome references.
    """
    from pycircuit.circuit.analysis import remove_row_col
    cir, pss, _pac, _rs = _loss_osc('series', Q=8.0, npts=480, a=0.25)
    m = cir.n - 1
    irn = pss.irefnode
    _v, info = pss.ppv()
    ## `samples_eq` IS `v_1` -- no reconstruction needed any more
    S = np.asarray(info['samples_eq'])[:, :m]
    Sv = np.asarray(info['samples'])[:, :m]
    X = np.asarray(pss.waveform[1], dtype=float)

    def mats(xf):
        Cr, Gr = remove_row_col((np.asarray(cir.C(xf), dtype=float),
                                 np.asarray(cir.G(xf), dtype=float)),
                                irn, pss.toolkit)
        return np.asarray(Cr, dtype=float), np.asarray(Gr, dtype=float)

    Cr0, _G0 = mats(X[:, 0])
    rows = [i for i in range(m) if not np.any(Cr0[i, :])]
    cols = [j for j in range(m) if not np.any(Cr0[:, j])]
    diff = [i for i in range(m) if i not in rows]
    assert rows and cols, 'this fixture must have an algebraic row to test'
    ## the branch row's C entry is negative -- the whole point of the sign
    assert Cr0[diff[-1], diff[-1]] < 0.0
    ## and `samples` must still be the STATE vector, structurally zero there
    assert np.max(np.abs(Sv[:, rows])) == 0.0, \
        'samples must remain C^T v_1, which annihilates the algebraic columns'

    scale = float(np.max(np.abs(S)))
    worst = worst_flipped = 0.0
    for k in range(0, S.shape[0], max(1, S.shape[0] // 8)):
        _Ck, Gk = mats(X[:, k])
        v1 = S[k]
        flipped = v1.copy()
        for i in rows:
            flipped[i] = -v1[i]
        worst = max(worst, float(np.linalg.norm(Gk[:, cols].T @ v1)))
        worst_flipped = max(worst_flipped,
                            float(np.linalg.norm(Gk[:, cols].T @ flipped)))
    assert worst / scale < 1e-12, \
        "Demir (24)'s algebraic rows are a pointwise constraint on v_1; the " \
        'fill violates it by %.3e (scaled)' % (worst / scale)
    ## and it discriminates -- the opposite sign is not merely worse, it is O(1)
    assert worst_flipped / scale > 1e-3, \
        'the constraint must REJECT the opposite sign, or it identifies ' \
        'nothing; it gave %.3e' % (worst_flipped / scale)


def _osc_for_deflation(Q=15.92, npts=400):
    """A van der Pol at a known `Q`, converged, for the deflated solve."""
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                  maxiterations=200)
    assert pss.converged
    return cir, pss


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_the_deflated_solve_is_capped_by_the_TANGENT_not_by_the_PPV():
    """⚠⚠ THE BORDER ROW'S ACCURACY DOES NOT ENTER THE ANSWER; the border
    COLUMN'S enters linearly -- BELOW THE REFINEMENT GATE.  Above it, since
    2026-09-22, neither does.

    The natural reading of a bordered solve is that it "consumes the null
    vectors", so its accuracy is capped by how well they are known.  That
    reading is HALF WRONG, and the half matters:
      - `v` (the PPV) enters only as the constraint row `v^T w = 0`.  For
        `alpha != 1` the system `(I - alpha M) y = b` is NONSINGULAR, so
        `y` is already determined by `b` alone; the border merely picks a
        well-conditioned route to it.  Any `v` not orthogonal to the null
        direction gives the SAME `y`.
      - `u` (the orbit tangent) enters the RECONSTRUCTION,
        `y = w + s u / (1 - alpha)`.  An error there is an error in the
        answer, and it passes straight through -- UNLESS the answer is
        then refined on the plain operator, which `_deflated_solve` does
        wherever `|1 - alpha| >= DEFLATION_REFINE_MIN` (1e-8): the
        refined answer is the operator's own and neither border vector
        is load-bearing any more.
    Measured through the shipped `_deflated_solve` (relative move in `y`):
        offset 1e-6 f0 (refined, |1 - alpha| 6e-6):
          perturbation   border ROW (v)    border COLUMN (u)
            1e-08          6.3e-11             3.7e-10
            1e-04          2.8e-10             2.5e-10
            1e-02          1.4e-10             5.0e-10
        offset 1e-10 f0 (below the gate, |1 - alpha| 6e-10):
            1e-08          1.2e-14             1.0e-08
            1e-04          5.9e-15             2.1e-04
            1e-02          9.8e-15             6.0e-03
    ⚠ SO A MORE ACCURATE PPV BUYS NOTHING HERE, and a more accurate orbit
    tangent buys everything only within 1e-8 of a harmonic.  Anyone
    tempted to tighten `ppv()`'s tolerance to improve a PAC result is
    optimising the wrong vector; and on a STAGED solve, whose discrete
    unit multiplier is displaced by O(h), the refinement is what makes
    the forward and adjoint solves agree (the staged-oscillator PAC test).
    """
    _cir, pss = _osc_for_deflation()
    fp = pss.factored_period()
    n = fp.width
    pac = PAC(_cir)
    rng = np.random.default_rng(1)
    b = rng.standard_normal(n).astype(complex)
    true_v, true_info = pss.ppv()
    true_v = np.asarray(true_v, float)
    true_u = np.asarray(true_info['tangent_pair'], float)
    orig = pss.ppv
    try:
        ## one part in 1e6 off the carrier (refined on the plain operator)
        ## and one part in 1e10 (below the refinement gate: the recovery alone)
        for offset, refined in ((1e-6, True), (1e-10, False)):
            alpha = np.exp(-2j * np.pi * (1.0 + offset))
            ref = pac._deflated_solve(pss, alpha, b)
            scale = float(np.linalg.norm(ref))
            assert scale > 1.0, \
                'the deflated answer is ~zero (%.3e), so the comparisons ' \
                'below would be vacuous' % scale
            for eps in (1e-8, 1e-4, 1e-2):
                d1 = rng.standard_normal(n)
                d1 /= np.linalg.norm(d1)
                d2 = rng.standard_normal(n)
                d2 /= np.linalg.norm(d2)
                vp = true_v + eps * np.linalg.norm(true_v) * d1
                up = true_u + eps * np.linalg.norm(true_u) * d2
                pss.ppv = lambda *a, **k: (vp, dict(true_info, tangent_pair=true_u))
                ev = float(np.linalg.norm(pac._deflated_solve(pss, alpha, b) - ref)) / scale
                pss.ppv = lambda *a, **k: (true_v, dict(true_info, tangent_pair=up))
                eu = float(np.linalg.norm(pac._deflated_solve(pss, alpha, b) - ref)) / scale
                pss.ppv = orig
                if refined:
                    assert ev < 1e-8 and eu < 1e-8, (offset, eps, ev, eu)
                    continue
                assert ev < 1e-11, \
                    'the border ROW now changes the answer (%.3e at eps=%.0e). ' \
                    'If that is real, the deflated solve has stopped being a ' \
                    'reformulation of a nonsingular system and the PPV\'s ' \
                    'accuracy has become load-bearing' % (ev, eps)
                assert eu > 0.05 * eps, \
                    'the border COLUMN no longer propagates linearly (%.3e at ' \
                    'eps=%.0e); the tangent is supposed to enter the ' \
                    'reconstruction directly' % (eu, eps)
                assert eu > 1e3 * max(ev, 1e-16), \
                    'the two vectors now matter comparably (row %.3e against ' \
                    'column %.3e at eps=%.0e); the asymmetry this test exists ' \
                    'to record is gone' % (ev, eu, eps)
    finally:
        pss.ppv = orig


def _osc_with_ladder(Q, nladder, nslow, npts=200):
    """A van der Pol at a target `Q` with an RC ladder whose first `nslow`
    sections have time constants STRADDLING the period.

    ⚠ THE STRADDLE IS THE WHOLE FIXTURE. A ladder whose sections all decay
    inside one step adds states without adding modes near the unit circle:
    it raises `m` and leaves the spectrum one cluster. That is what makes
    `nslow` and `m` separable here, and separating them is the point --
    the first version of this measurement confounded them and reported a
    flat iteration count at every `(Q, m)`.
    """
    circuit.default_toolkit = circuit.numeric
    tper = 2.0 * np.pi
    mu = 1.0 / (2.0 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    prev = 'v'
    for j in range(nladder):
        nd = 'p%d' % j
        cir.add_node(nd)
        tau = (tper * 10.0 ** (-1.0 + 2.0 * j / max(nslow - 1, 1))
               if j < nslow else tper * 1e-4)
        cir['r%d' % j] = R(prev, nd, r=1e3)
        cir['c%d' % j] = C(nd, gnd, c=tau / 1e3)
        prev = nd
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-11)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=200)
    assert pss.converged, 'Q=%g/nslow=%d did not converge' % (Q, nslow)
    return cir, pss


def _bordered_gmres_iterations(pss):
    """GMRES iterations for the bordered `(I - M) w = b` on an oscillator."""
    import scipy.sparse.linalg as spla
    fp = pss.factored_period()
    n = fp.width
    with quiet(AccuracyWarning):
        v_, info = pss.ppv()
    v = np.asarray(v_, dtype=float)
    u = np.asarray(info['tangent_pair'], dtype=float)

    def mv(z):
        z = np.asarray(z)
        w, s = z[:n], z[n]
        return np.concatenate((w - np.asarray(fp.matvec(w)) + s * u,
                               [float(v @ w)]))

    A = spla.LinearOperator((n + 1, n + 1), matvec=mv, dtype=float)
    rng = np.random.default_rng(0)
    b = np.concatenate((rng.standard_normal(n), [0.0]))
    its = [0]
    spla.gmres(A, b, rtol=1e-10, restart=min(n + 1, 200), maxiter=50,
               callback=lambda *a: its.__setitem__(0, its[0] + 1),
               callback_type='pr_norm')
    return its[0]


def test_krylov_cost_ignores_Q_and_tracks_the_SLOW_NODE_COUNT():
    """⚠⚠ THE HIGH-Q WORRY IS FALSIFIED, AND THE REAL DRIVER IS SOMETHING
    ELSE.

    The open question was whether a matrix-free shooting solve collapses
    at high `Q`: the multipliers crowd the unit circle, so `I - M` has its
    spectrum crowding zero, and GMRES was expected to need `O(m)`
    iterations exactly where large `m` makes matrix-free worth having.

    Measured at FIXED `m = 32`, sweeping the number of ladder sections
    whose time constant straddles the period:

        nslow      0    4    8   14   22   30
        Q =   8    4    8   11   16   23   29
        Q = 256    5    8   11   17   23   30

    ⚠ A 32x CHANGE IN `Q` MOVES THE COUNT BY AT MOST ONE. Iterations
    track `nslow` -- roughly `1 + nslow` -- and ignore both `Q` and `m`.
    Krylov iteration count is set by the number of DISTINCT eigenvalue
    clusters, which is a spectral-spread property, not by conditioning,
    which is what `Q` controls.

    ⚠ SO THE OPERATIONAL RULE INVERTS: matrix-free is SAFE on a high-Q
    oscillator and degrades on a circuit with many SLOW NODES, at any `Q`.
    A designer's high-Q tank costs nothing here; a bias network with a
    dozen long time constants costs linearly.

    ⚠ AND `|lambda| > 0.9` IS A POOR PROXY for the driver -- it counted
    1/2/2/3/4/5 across that sweep while iterations went 4/8/11/16/23/29.
    The count of near-unit multipliers above a threshold is not the same
    as the number of distinct clusters, and only the latter predicts.
    """
    ## m = 16 and a coarser grid than the sweep above: the CONTRAST is
    ## what is being pinned, not the absolute counts
    its = {}
    for Q in (8.0, 256.0):
        for nslow in (0, 14):
            _cir, pss = _osc_with_ladder(Q, 14, nslow)
            its[(Q, nslow)] = _bordered_gmres_iterations(pss)

    ## 1. slow nodes cost, and cost a lot
    for Q in (8.0, 256.0):
        assert its[(Q, 14)] > 2 * its[(Q, 0)], \
            'Q=%g: a fully slow ladder (%d iterations) no longer costs ' \
            'materially more than a fully fast one (%d) at the same m -- ' \
            'the fixture has stopped separating the two' \
            % (Q, its[(Q, 14)], its[(Q, 0)])

    ## 2. ⚠ AND Q DOES NOT. This is the falsification, and it is the
    ## assertion that would break if the high-Q worry were real.
    for nslow in (0, 14):
        d = abs(its[(8.0, nslow)] - its[(256.0, nslow)])
        assert d <= 3, \
            'at nslow=%d the iteration count moved by %d across a 32x ' \
            'change in Q (%d against %d). Krylov cost is supposed to be ' \
            'insensitive to Q; if that has changed, the matrix-free route ' \
            'is no longer safe on high-Q oscillators and the roadmap\'s ' \
            'conclusion needs re-measuring' \
            % (nslow, d, its[(8.0, nslow)], its[(256.0, nslow)])


def test_floquet_modes_are_genuinely_periodic():
    """A9's prerequisite: the Floquet pairs, with the periodic part.

    ⚠⚠ **`|λ₂|` ALONE IS NOT ENOUGH FOR THE ORBITAL SPECTRUM, BY THE
    SOURCE'S OWN STATEMENT.** Traversa & Bonani (TCAS-I 2011) make `S_yy`
    a sum of Lorentzians weighted by `C_lhj` (their eq 22), which is built
    from the **Fourier coefficients of `u_l(t)` and `v_l(t)ᵀB(t)`** — and
    their §III says in terms that *"a major role in the C and D
    coefficients is also played by the Floquet eigenvectors, which could
    determine large orbital fluctuations contributions even when the
    Floquet exponents are not near to zero."* So the exponents do not
    order the result and the eigenvectors are not optional.

    ⚠ THE GATE IS FLOQUET'S THEOREM ITSELF, which needs no reference:
    the solution is `p_l(t)·exp(μ_l t)` with `p_l` **T-periodic**, so
    `p_l(T) = p_l(0)`. That holds only if `λ_l`, `μ_l = log(λ_l)/T` and
    the propagation are all consistent — a wrong multiplier breaks it
    even when the eigenvector residual is clean.

    Measured on van der Pol: periodicity 2.8e-15 (the unit mode) and
    5.9e-15 (the amplitude mode), with eigenvector residuals 9.2e-16 and
    4.0e-16.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged

    ## ⚠ THE DEFAULT PATH FIRST. `nmodes=None` returns every non-null mode
    ## and is the documented default; it raised `int(None)` for an hour
    ## because this test only ever passed a number.
    modes = pss.floquet_modes()
    assert len(modes) == 2, \
        'the default (all non-null modes) returned %d, expected 2' % len(modes)
    modes = pss.floquet_modes(nmodes=2)
    assert len(modes) == 2, 'expected two non-null modes, got %d' % len(modes)

    ## the phase mode is the unit multiplier, and it must come first
    assert abs(abs(modes[0]['lam']) - 1.0) < 1e-9, \
        'the leading multiplier is %.12f, not 1 — an autonomous ' \
        'oscillator must carry the phase mode' % abs(modes[0]['lam'])
    ## and the second is the amplitude mode, strictly inside
    assert abs(modes[1]['lam']) < 1.0 - 1e-6, \
        'the second multiplier is not inside the unit circle (%.12f)' \
        % abs(modes[1]['lam'])

    for k, md in enumerate(modes):
        assert md['residual'] < 1e-10, \
            'mode %d eigenvector residual %.3e' % (k, md['residual'])
        P = md['p']
        per = float(np.linalg.norm(P[:, -1] - P[:, 0])) \
            / max(float(np.linalg.norm(P[:, 0])), 1e-300)
        assert per < 1e-9, \
            'mode %d: p(T) differs from p(0) by %.3e, so the propagated ' \
            'vector is NOT Floquet-periodic. Either the multiplier, the ' \
            'exponent log(lam)/T, or the propagation disagrees with the ' \
            'other two — this is the check that catches a wrong lambda ' \
            'even when the eigenvector residual is clean' % (k, per)

    ## ⚠⚠ AND THE STATE-BLOCK PAIR MUST BE BIORTHONORMAL, which the
    ## periodicity check above CANNOT see (it is scale-free). Under gear
    ## the width-n normalisation left q(0)^T p(0) = 1.324 on the width-m
    ## block, and the orbital covariance built from these parts was too
    ## large by exactly 1.324^2 against two independent routes.
    ##
    ## ⚠⚠⚠ THE INNER PRODUCT IS `C`-WEIGHTED, AND THIS ASSERTION USED THE
    ## UNWEIGHTED ONE UNTIL 2026-09-07.  The variational DAE conserves
    ## `q(t)^T C(t) p(t)`, not `q(t)^T p(t)`; the two are the SAME NUMBER on
    ## a unit-reactance fixture, which is every fixture this file had, so the
    ## error was invisible here and cost a factor of `C^2` in every orbital
    ## covariance (measured 16.04x too large at `c = 4` against the
    ## independent Lyapunov reference, and right only at `c = 1`).
    ## ⚠ On van der Pol the reduced `C` carries a NEGATIVE entry for the
    ## inductor, so the correctly-normalised modes now read
    ## `q^T p = -1.001` — the sign is the tell that the old assertion was
    ## measuring a different bilinear form, not a scale.
    ## The `c != 1` sweep lives in
    ## `test_the_orbital_spectrum_is_a_lorentzian_of_half_width_f_amp`;
    ## this one pins the invariant itself, AROUND THE CYCLE, because
    ## conservation is the actual claim and `t = 0` alone would not test it.
    ##
    ## ⚠⚠⚠ AND `q` IS STORED IN REVERSE TIME ORDER RELATIVE TO `p`, WHICH
    ## WAS UNDOCUMENTED AND IS MEASURED HERE.  The adjoint comes from a
    ## BACKWARD replay.  Pairing at the same index gives nonsense --
    ## `q[j]^T C p[j]` = 1.000, -0.022, -1.000, 1.000, -0.999 across the
    ## cycle -- while `q[N-1-j]^T C p[j]` = 1.000, 0.990, 0.999, 1.000,
    ## 1.000 is the conserved invariant.  The two agree only at `j = 0`,
    ## `N/2` and `N-1` by the orbit's symmetry, and index 0 is the ONLY
    ## place any shipped code looked, which is why nothing caught it.
    ##
    ## ⚠ FLAGGED, NOT RESOLVED: `orbital_correlation` matches an INDEPENDENT
    ## Lyapunov reference at `c = 0.25/1/4` to 0.4 %, so its Fourier path is
    ## self-consistent with this storage.  But any consumer that pairs
    ## `q[:, k]` with `p[:, k]` AT THE SAME k is wrong -- including the
    ## definition-route integral inside
    ## `test_orbital_correlation_is_gated_three_ways`, which may be masked
    ## by van der Pol's half-wave symmetry averaging the mismatch out over
    ## the cycle.  Whether the reversal is intentional storage or a latent
    ## defect is OPEN; it needs an ASYMMETRIC orbit to separate.
    x0r = np.delete(np.asarray(pss.waveform[1], dtype=float)[:, 0],
                    pss.irefnode)
    Cm = np.asarray(pss._C_at(x0r), dtype=float)
    for k, md in enumerate(modes):
        c = complex(np.vdot(md['q'][:, 0], Cm @ md['p'][:, 0]))
        assert abs(c - 1.0) < 1e-9, \
            'mode %d: q(0)^T C p(0) = %.6f%+.6fj on the state block, not 1 -- ' \
            'any covariance assembled from these parts is off by |c|^2' \
            % (k, c.real, c.imag)
        ## ⚠⚠ AROUND THE CYCLE, AT THE SAME INDEX -- RESTORED 2026-09-07 after
        ## being removed the same day.  This block once asserted
        ## `q[j]^T p[j] = 1`, which was wrong twice (no `C` weighting; and it
        ## read a spread of 2.0).  A reversed pairing `q[N-1-j]` then gave
        ## ~1e-2 that did NOT converge with refinement, and was recorded as
        ## "q is stored in reverse time order, correspondence unknown".
        ##
        ## ⚠ BOTH READINGS WERE ARTIFACTS OF ONE DEFECT: the replayed adjoint
        ## was `C^T q`, not `q` (see the `C^-T` block in `floquet_modes`).
        ## With `C = diag(1, -1)` that flips one component, which on a
        ## half-wave symmetric orbit is exactly the relation between `q(t)`
        ## and `q(T - t)` -- a sign flip read as a time reversal.  With the
        ## right vector the invariant holds at the SAME index (measured
        ## spread 4.2e-04 on this fixture) and NOT at the reversed one (2.0).
        ## ⚠ `C(t_j)` PER SAMPLE, not `C(0)` for every sample.  Demir (DAEs and
        ## Colored Noise Sources, Remark 2.1) states the convention as
        ## `v_j^T(t) C(t) u_i(t) = delta` and derives the invariant for the
        ## CHARGE-based forward `d/dt(C x) = -G x`; with `C` time-varying
        ## (a varactor, a forward-biased junction) `C(0)` here would falsely
        ## fail on an orbit the code handles correctly.  On these linear-
        ## reactance fixtures the two coincide, which is why it passed either
        ## way -- the same blindness as everything else today.
        Pm_, Qm_ = md['p'], md['q']
        _nn = min(Pm_.shape[1], Qm_.shape[1])
        _Wc = np.delete(np.asarray(pss.waveform[1], dtype=float),
                        pss.irefnode, axis=0)
        cyc = [abs(complex(np.vdot(
                   Qm_[:, j],
                   np.asarray(pss._C_at(_Wc[:, min(j, _Wc.shape[1] - 1)]),
                              dtype=float) @ Pm_[:, j])) - 1.0)
               for j in range(0, _nn - 1, max(1, _nn // 8))]
        assert max(cyc) < 5e-3, \
            'mode %d: q(t)^T C p(t) drifts from 1 around the cycle by %.3e ' \
            'at the SAME index. If this has regressed to ~2, the replayed ' \
            'adjoint is being used as q without the C^-T transform again' \
            % (k, max(cyc))

    ## ⚠ AND THE NULL MODES MUST BE ABSENT. A DAE monodromy has exact
    ## zeros; asked for more modes than exist, it must not pad with them.
    many = pss.floquet_modes(nmodes=10)
    assert all(abs(md['lam']) > 1e-12 for md in many), \
        'a null (annihilated algebraic) multiplier was returned as a mode'


def test_the_orbit_is_read_full_width_past_the_reference_node():
    """⚠ `PSS.waveform` is FULL width -- the reference row is in it -- and
    two readers took it for the reduced state and inserted the reference's
    zero a second time, reading every unknown past the reference one slot
    late (found 2026-09-26 when a DC source entered an oscillator):

      * `PAC._phase_mode_split` computed the orbit tangent from it: on a
        van der Pol with an idle 3 V source the phase mode's alignment read
        0.57 and `modal_spectrum` refused the circuit.  Plain van der Pol
        hides it -- the inductor current it misread is ~0 at the phase
        anchor.
      * `PAC._cy_cycle_averaged` (`pnoise(modulated=True)`) evaluated
        `CY` at the shifted states: a source controlled by V(v, b) read
        `b` as 0.
    Both now read `_orbit_states`, right whichever width it is given."""
    _c, pss, pac, ov = _orbit_modulated_vdp('lorentz')
    with quiet():
        modes = pss.floquet_modes()
    k, _orb = pac._phase_mode_split(pss, modes, 'test')
    assert abs(abs(complex(modes[k]['lam'])) - 1.0) < 1e-9
    w = 2.0 * np.pi / float(pss.period)
    fp = pss.factored_period()
    Cs = pac._noise_components(pss).cy_at_states(w)
    hs = pac._period_weights(np.asarray(fp.times, dtype=float), Cs.shape[0],
                             float(fp.T), pss)
    ref = np.tensordot(hs, Cs, axes=1) / float(hs.sum())
    got = pac._cy_cycle_averaged(pss, w)
    assert np.max(np.abs(got - ref)) <= 1e-12 * np.max(np.abs(ref)), \
        (got, ref)


def test_the_raw_pair_dc_is_the_consistent_dc_times_1p5_s():
    """✅ THE LAST RESIDUE OF §0l, CLOSED: the raw pair block's DC content is
    the consistent object's times `1.5 s`, with `s` the pair-consistency
    scale read off the stored second blocks (`samples_pair[:, m:] = w2`,
    `samples[:, m:] = w2 / s`).  `s = 2/3` exactly for an isochronous pair
    (`w2 = -w1/3`), so the raw block's DC error `1.5 s - 1` is +8.1% on the
    bias-sensitive core (`s = 0.7207`) and +0.014% on the divider -- and
    the divider's node-v mean is 4e-6 |v|, so 0.014% of it is the 4e-11
    absolute that had been recorded as an unexplained exactness.  A units
    mix (absolute on one side, relative on the other), the review
    session's shape 0i.  Pinned on the bias core's inductor row and the
    series-loss tank's two rows: `mean(raw)/mean(consistent) = 1.5 s` to
    2e-4.
    """
    circuit.default_toolkit = circuit.numeric

    def bias():
        mu = 1.0 / (2.0 * np.pi * 8.0)
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0)
                           + 0.3 * u ** 2)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with quiet(AccuracyWarning):
            pss.solve(period=6.731, timestep=6.731 / 400,
                      x0=np.array([2.0, 0.0]), maxiterations=300)
        return cir, pss, [1]

    def lossy():
        cir = SubCircuit()
        cir.add_node('v')
        cir.add_node('x')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', 'x', L=1.0)
        cir['Rs'] = R('x', gnd, r=0.2)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: 1.0 * (u - u ** 3 / 3.0)
                           + 0.25 * (u ** 2 - 2.0))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        pss = PSS(cir, method='gear', reltol=1e-12)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        with quiet(AccuracyWarning):
            pss.solve(period=6.66, timestep=6.66 / 480, x0=x0,
                      maxiterations=200)
        return cir, pss, [0, 2]
    expect_s = {'bias': 0.7207, 'lossy': 0.6656}
    for name, build in (('bias', bias), ('lossy', lossy)):
        cir, pss, rows = build()
        assert pss.converged
        m = cir.n - 1
        _v, info = pss.ppv()
        h = np.diff(np.asarray(info['times'], dtype=float))
        n = len(h)
        T = float(pss.period)
        Sp = np.asarray(info['samples_pair'], dtype=float)[:n]
        Sc = np.asarray(info['samples'], dtype=float)[:n]
        w2r, w2c = Sp[:, m:], Sc[:, m:]
        mask = np.abs(w2c) > 1e-3 * np.abs(w2c).max()
        s = float(np.median(w2r[mask] / w2c[mask]))
        assert abs(s - expect_s[name]) < 2e-3, \
            '%s: pair-consistency scale s = %.4f, expected %.4f' \
            % (name, s, expect_s[name])
        mr = (Sp[:, :m] * h[:, None]).sum(0) / T
        mc = (Sc[:, :m] * h[:, None]).sum(0) / T
        for j in rows:
            ratio = mr[j] / mc[j]
            assert abs(ratio / (1.5 * s) - 1.0) < 2e-4, \
                '%s row %d: mean(raw)/mean(consistent) = %.6f against ' \
                '1.5 s = %.6f' % (name, j, ratio, 1.5 * s)


def test_the_ppv_waveform_matches_a_pulse_isf_over_the_whole_period():
    """The PPV checked as a WAVEFORM, by kicks at phases around the orbit.

    ⚠ THIS EXISTS TO RETIRE A WEAKNESS THE t=0 GATE NAMES ABOUT ITSELF.
    `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has` kicks
    in a RANDOM direction, and its own docstring says why that is fragile:
    "a random direction in TWO dimensions is ~71% tangential, so the phase
    signal dominates ... THAT PROTECTION SCALES AS 1/sqrt(m) AND VANISHES ON
    A REAL CIRCUIT ... This gate is sound at m = 2 and would not be at
    m = 20, with nothing in it changing."

    This gate kicks along COORDINATE directions at ten phases spread over
    the period.  There is no random direction in it, so it carries no
    `1/sqrt(m)` dependence, and it exercises `info['samples']` -- the PPV
    over the orbit -- rather than the single vector at `t = 0`.

    Measured (van der Pol, mu = 1, 400 points), 20 independent pulse
    experiments, worst |1 - measured/predicted| = **4.2e-03**::

        t/T     e0 measured     e0 predicted    ratio
        0.000   +8.113272e-02   +8.145265e-02   0.9961
        0.201   -7.161729e-01   -7.162341e-01   0.9999
        0.501   -7.836864e-02   -7.868975e-02   0.9959
        0.702   +7.144444e-01   +7.144991e-01   0.9999

    ⚠⚠ THE HALF-WAVE ANTISYMMETRY IS THE SELF-CHECK, AND NOTHING IN THE
    MEASUREMENT IMPOSES IT.  Van der Pol is half-wave symmetric, so its ISF
    inherits `Gamma(t + T/2) = -Gamma(t)`.  The pulse experiments at `t` and
    at `t + T/2` are entirely independent transients -- different initial
    states, different trajectories -- so agreement between them is evidence
    the harness is sound, not an identity it was built to satisfy.

    ⚠ THE INDEX CONVENTION IS PINNED, NOT ASSUMED.  `info['samples']` comes
    from a REVERSE replay, so whether `samples[k]` is `t_k` or `t_{N-1-k}`
    is exactly the off-by-one that has bitten this arc before (the
    sideband-fold abscissa, the conjugation).  It is settled here by the
    normalisation `v(t).xdot(t) = 1`, which holds at every k for the
    forward reading and gives 0.42 / -1.04 for the reversed one -- so the
    assertion below would FAIL on an index flip rather than absorb it.

    An independent session implementing Levantino's reference pulse method
    reported the same waveform (+8.11e-2, -3.84e-1, -7.17e-1, -3.36e-1,
    -1.56e-1 at t/T = 0 .. 0.4); this reproduces those numbers from a
    separately written harness.
    """
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric

    npts = 400
    cir, pss, v, info = _vdp_ppv(npts)
    m = cir.n - 1
    irn = pss.irefnode
    T = pss.period
    Xf = np.asarray(pss.waveform[1], dtype=float)
    Xr = np.delete(Xf, irn, axis=0)
    S = np.asarray(info['samples'], dtype=float)
    ts = np.asarray(info['times'], dtype=float)
    nint = Xr.shape[1] - 1                     # intervals in the period

    def xdot_at(k):
        h = T / nint
        return (Xr[:, (k + 1) % nint] - Xr[:, (k - 1) % nint]) / (2.0 * h)

    ## (1) PIN THE ORDERING.  `v(t).xdot(t) = 1` is the normalisation, which
    ## makes it the right instrument for an INDEX question and the wrong one
    ## for a correctness question -- it is used only for the former.
    probe = (0, 50, 100, 200, 300)
    fwd = [float(S[k][:m] @ xdot_at(k)) for k in probe]
    rev = [float(S[len(S) - 1 - k][:m] @ xdot_at(k)) for k in probe]
    assert max(abs(z - 1.0) for z in fwd) < 5e-3, \
        'samples[k] <-> t_k should satisfy v.xdot = 1, got %r' % (fwd,)
    assert max(abs(z - 1.0) for z in rev) > 0.1, \
        'the REVERSED reading also satisfies the normalisation (%r), so this ' \
        'gate cannot tell an index flip from the truth -- it must' % (rev,)

    def integrate(xi, ppp=2000):
        tran = Transient(cir, toolkit=circuit.numeric, reltol=1e-9,
                         iabstol=1e-13, vabstol=1e-11)
        with quiet():
            res = tran.solve(refnode=gnd, tend=T, timestep=T / ppp, x0=xi)
        return np.asarray(res.x, dtype=float)[:, -1]

    ## (2) THE WAVEFORM.  Phases chosen in half-period PAIRS so the same
    ## transients serve the antisymmetry check below -- no extra cost.
    eps = 1e-5
    half = nint // 2
    ks = [0, 40, 80, half, half + 40, half + 80]
    gamma = {}
    worst = 0.0
    for k in ks:
        ref = integrate(Xf[:, k].copy())
        xd = xdot_at(k)
        for j in range(m):
            d = np.zeros(m)
            d[j] = 1.0
            dr = np.concatenate((d[:irn], np.zeros(1), d[irn:]))
            dx = np.delete(integrate(Xf[:, k].copy() + eps * dr) - ref, irn)
            meas = float(dx @ xd) / float(xd @ xd) / eps
            pred = float(S[k][:m] @ d)
            gamma[(k, j)] = meas
            assert abs(pred) > 1e-3, \
                'the PPV is ~0 at t/T=%.3f along e%d, so this point is a ' \
                'zero-vs-zero pass' % (ts[k] / T, j)
            worst = max(worst, abs(1.0 - meas / pred))
    assert worst < 1.5e-2, \
        'the PPV waveform disagrees with the pulse ISF by %.3e at worst' % worst

    ## (3) HALF-WAVE ANTISYMMETRY of the MEASURED waveform -- independent
    ## transients, so this is evidence about the harness, not an identity.
    ## ⚠ NORMALISED BY THE WAVEFORM'S PEAK, NOT BY THE LOCAL VALUE.  The
    ## first version of this divided by `max(|a|,|b|)` and reported a 10%
    ## violation -- all of it from the `e1` pair near a ZERO CROSSING
    ## (+5.29e-2 against -4.69e-2), where a small denominator inflates a
    ## small absolute difference.  That is a defect in the measure, not in
    ## the waveform: a relative error against a quantity passing through
    ## zero is not a statement about agreement.
    scale = max(abs(z) for z in gamma.values())
    anti = 0.0
    for k in (0, 40, 80):
        for j in range(m):
            a, b = gamma[(k, j)], gamma[(k + half, j)]
            anti = max(anti, abs(a + b) / scale)
    assert anti < 2e-2, \
        'van der Pol is half-wave symmetric so its ISF must obey ' \
        'Gamma(t+T/2) = -Gamma(t); worst violation %.3e of the peak' % anti


def _dense_lam2_of(fp, deflate=1e-6):
    """`lam2` from the DENSE spectrum of the same operator -- `n` matvecs and
    `eigvals`, the identical route `PSS.ppv` already takes for `dirk`/`full`,
    with the identical selection rule.  It cannot be influenced by the Arnoldi
    it is the reference for."""
    n = fp.width
    M = np.column_stack([np.asarray(fp.matvec(e), float).ravel()
                         for e in np.eye(n)])
    lams = np.linalg.eigvals(M)
    keep = np.real(lams)[np.abs(lams - 1.0) > deflate]
    return (float(max(np.max(keep), 0.0)) if keep.size else 0.0), lams


def _arnoldi_lam2_and_ritz_residual(fp, kk, deflate=1e-6):
    """`PSS.ppv`'s Arnoldi, replicated exactly (same seed, same selection),
    plus the per-pair Ritz residual `|h_{k+1,k}| |y_i[last]|` for the pair it
    selects -- which needs no extra matvec and is not currently computed."""
    n = fp.width
    kk = int(min(n, kk))
    rng = np.random.default_rng(12345)
    q0 = rng.standard_normal(n)
    q0 = q0 / np.linalg.norm(q0)
    Qb, H = [q0], np.zeros((kk + 1, kk))
    for j in range(kk):
        wj = Qb[j] - np.asarray(fp.matvec(Qb[j]))
        for i in range(j + 1):
            H[i, j] = float(Qb[i] @ wj)
            wj = wj - H[i, j] * Qb[i]
        H[j + 1, j] = float(np.linalg.norm(wj))
        if H[j + 1, j] < 1e-13:
            kk = j + 1
            break
        Qb.append(wj / H[j + 1, j])
    theta, Y = np.linalg.eig(H[:kk, :kk])
    lams = 1.0 - theta
    mask = np.abs(lams - 1.0) > deflate
    if not mask.any():
        return 0.0, float('nan')
    lam2 = float(max(np.max(np.real(lams)[mask]), 0.0))
    idx = int(np.where(mask)[0][int(np.argmax(np.real(lams)[mask]))])
    res = abs(H[kk, kk - 1]) * abs(Y[kk - 1, idx]) / max(
        float(np.linalg.norm(Y[:, idx])), 1e-300)
    return lam2, float(res)


@pytest.mark.slow
def test_the_ppv_takes_the_dense_spectrum_when_it_can_afford_it():
    """`lam2` comes from the SPECTRUM below `FLOQUET_DENSE_LIMIT`, not an Arnoldi.

    ⚠⚠ THIS TEST BEGAN LIFE ASSERTING THE DEFECT.  It was written to pin a
    KNOWN GAP -- `PSS.ppv` reported `info['second_multiplier']` and `info['Q']`
    from a `k = PPV_RITZ_BASIS = 12` Arnoldi on `I - M`, and got them wrong by
    4x-19x once slow nodes crowded the unit root.  It carried a message telling
    whoever fixed it to delete the "wrong" branch, and that is what happened.
    The measurement it was built on is below, because the fix is only as good
    as the reason for it.

    Against the dense spectrum of the SAME operator, scored in the GAP because
    `Q ~ 1/(1 - lam2)`, on `_osc_with_ladder(16, 14, nslow)`::

        nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense / Arnoldi
         <=11   0.995706203    0.995706197      1.000        232 / 232
           12   0.996324417    1.000114048     -0.031        271 / inf
           13   0.996818781    0.942674586     18.020        313 / 16.9
           14   0.997220139    0.999318472      0.245        359 / 1467

    It overturned two claims `ppv` recorded as justification for the cap:
    that a truncated `lam2` is a CAUCHY LOWER BOUND so the near-unit warning
    can only under-fire (true for a NORMAL `M`; a circuit monodromy is not one,
    and the error above is **not one-signed**), and that the selection failure
    was NOT LIVE on a circuit monodromy (measured on the eigenvector-
    conditioning axis; the trigger is the CLUSTER COUNT, a different axis, which
    `_osc_with_ladder` varies by construction).

    ⚠ AND RAISING THE CONSTANT WAS MEASURED NOT TO BE THE FIX, which is why the
    fix is a route change and not a bigger number::

        ladder/nslow   n    k=12    k=16    k=20
           14 / 14     32   0.245   1.000   1.000
           20 / 20     44   0.277   0.279   1.000
           26 / 26     56   2.313   0.410   0.265

    The basis has to grow with the problem; a constant cannot. `k = 16` passes
    the fixture above and ships the same defect on a longer ladder.

    ⚠ WHAT IS STILL NOT FIXED, and the warning says so: above
    `FLOQUET_DENSE_LIMIT` the truncated Arnoldi is all there is, and that is
    where it is least trustworthy -- a big circuit is the one likely to carry
    many slow nodes. Lai's 64-gated-capacitor DCO is >500 equations. The
    remaining answer is a per-pair Ritz residual gate; the roadmap has it.
    """
    circuit.default_toolkit = circuit.numeric

    ## (1) THE REFERENCE MUST BE SOUND BEFORE AGREEMENT WITH IT MEANS ANYTHING.
    fps = {}
    for nslow in (11, 12, 13, 14):
        _cir, pss = _osc_with_ladder(16.0, 14, nslow)
        fp = pss.factored_period()
        ld, lams = _dense_lam2_of(fp)
        unit = float(np.min(np.abs(lams - 1.0)))
        assert unit < 1e-9, \
            'nslow=%d: the dense unit root sits at |lam-1| = %.2e, not ' \
            'decisively inside the 1e-6 deflation -- the REFERENCE is then ' \
            'as ambiguous as the thing it judges' % (nslow, unit)
        assert 0.99 < ld < 1.0, \
            'nslow=%d: dense lam2 = %.6f, outside the near-unit regime this ' \
            'test is about' % (nslow, ld)
        fps[nslow] = (fp, ld, pss)

    ## (2) AND THE FIXTURE MUST STILL CROWD THE UNIT ROOT, or there is no
    ## cluster count to have been the trigger.
    _ld14, lams14 = _dense_lam2_of(fps[14][0])
    above = int(np.sum(np.abs(lams14) > 0.9))
    assert above >= 4, \
        'nslow=14 has only %d multipliers above 0.9; the ladder has stopped ' \
        'manufacturing the crowding this test is about' % above

    ## (3) THE FIX: `ppv` agrees with the spectrum at EVERY nslow, including
    ## the three that used to be wrong, and says which route it took.
    for nslow in (11, 12, 13, 14):
        fp, ld, pss = fps[nslow]
        with quiet(AccuracyWarning):
            _v, info = pss.ppv()
        la = float(info['second_multiplier'])
        assert info['second_multiplier_route'] == 'dense', \
            'nslow=%d: n = %d is inside FLOQUET_DENSE_LIMIT = %d, so this ' \
            'must come from the spectrum, not route %r' \
            % (nslow, fp.width, PSS.FLOQUET_DENSE_LIMIT,
               info['second_multiplier_route'])
        ratio = (1.0 - la) / (1.0 - ld)
        assert abs(ratio - 1.0) < 1e-9, \
            'nslow=%d: gap ratio %.6f against the dense spectrum. The k=12 ' \
            'Arnoldi gave -0.031 / 18.020 / 0.245 at nslow 12/13/14; if this ' \
            'has come back, the dense route is no longer being taken.' \
            % (nslow, ratio)

    ## (4) ⚠ NEUTER = THE OLD ROUTE.  Forcing the truncated path back must
    ## reproduce the recorded failure, or this test no longer guards the
    ## reason the fix exists.  Both signs, and the `lam2 > 1` case.
    bad = {}
    for nslow in (12, 13, 14):
        fp, ld, pss = fps[nslow]
        pss.FLOQUET_DENSE_LIMIT = 4          # below n = 32, so Arnoldi again
        ## ⚠ AND THE BASIS MUST BE STARVED TOO, which is itself a result: the
        ## Ritz-residual gate GROWS `k` until the pair certifies, so the
        ## truncated path now gets these right on its own and the old defect
        ## is unreachable without disabling both mechanisms.  See
        ## `test_the_truncated_lam2_is_gated_on_its_own_ritz_residual`.
        pss.PPV_RITZ_MAX_BASIS = 12
        ## `ppv` holds no cache -- it recomputes -- so re-calling it under
        ## the lowered limit really does take the other branch, which the
        ## route assertion below checks rather than assumes.
        with quiet(AccuracyWarning):
            _v, info = pss.ppv()
        assert info['second_multiplier_route'] == 'arnoldi', \
            'the neuter did not reach the truncated path'
        la = float(info['second_multiplier'])
        bad[nslow] = (la, (1.0 - la) / (1.0 - ld))
    for nslow in (12, 13, 14):
        assert abs(bad[nslow][1] - 1.0) > 0.5, \
            'NEUTER: nslow=%d no longer fails on the truncated path (gap ' \
            'ratio %.3f), so the defect this fix removes has gone somewhere ' \
            'else and the fix is unguarded' % (nslow, bad[nslow][1])
    assert bad[13][1] > 1.0 and bad[14][1] < 1.0, \
        'the failure was NOT ONE-SIGNED -- an under-estimate at 13 (%.3f) and ' \
        'an over-estimate at 14 (%.3f). That is what killed the Cauchy ' \
        '"can only under-fire" claim, and it is the half worth keeping.' \
        % (bad[13][1], bad[14][1])
    assert bad[12][0] > 1.0, \
        'and nslow=12 returned lam2 = %.9f > 1 -- a spurious UNSTABLE ' \
        'multiplier, which `Q` reports as inf' % bad[12][0]

    ## (5) ⚠ AND THE CONSTANT WAS NOT THE FIX.  `k = 16` is exact on the
    ## fixture above and fails on a longer ladder, so a bigger
    ## `PPV_RITZ_BASIS` would have passed this test and shipped the defect.
    _c20, p20 = _osc_with_ladder(16.0, 20, 20)
    fp20 = p20.factored_period()
    ld20, _l20 = _dense_lam2_of(fp20)
    la20_16, res20_16 = _arnoldi_lam2_and_ritz_residual(fp20, 16)
    ratio20 = (1.0 - la20_16) / (1.0 - ld20)
    assert abs(ratio20 - 1.0) > 0.5, \
        'k=16 now gets the longer ladder right too (gap ratio %.3f). If that ' \
        'holds at nladder=26 as well, a constant really would have been ' \
        'enough and this warning can go; it did not when measured (0.410).' \
        % ratio20

    ## (6) THE DIAGNOSTIC THAT WOULD MAKE THE REMAINING TRUNCATED PATH SAFE,
    ## recorded because that path still exists above the limit: right answers
    ## and wrong ones separate by thirteen orders.
    res_ok = _arnoldi_lam2_and_ritz_residual(fps[11][0], 12)[1]
    res_bad = [_arnoldi_lam2_and_ritz_residual(fps[n][0], 12)[1]
               for n in (12, 13, 14)]
    assert res_ok < 1e-6 < min(res_bad), \
        'the per-pair Ritz residual no longer separates right (%.2e) from ' \
        'wrong (%r); then the roadmap\'s recommendation for n > ' \
        'FLOQUET_DENSE_LIMIT needs re-measuring' % (res_ok, res_bad)


@pytest.mark.slow
def _align_mode(ref, got):
    """`got`'s `p`, `q` rescaled onto `ref`'s: a mode is defined up to a scalar
    `a` (`p -> a p`, `q -> q / conj(a)` keeps `q^T C p = 1`)."""
    a = complex(np.vdot(got['p'][:, 0], ref['p'][:, 0])
                / np.vdot(got['p'][:, 0], got['p'][:, 0]))
    return got['p'] * a, got['q'] / np.conj(a)


def test_floquet_modes_above_the_dense_limit_are_the_dominant_ritz_certified_modes():
    """Above `FLOQUET_DENSE_LIMIT` (2026-09-25; refused outright before)
    `floquet_modes(nmodes=k)` returns the k DOMINANT modes from a
    Ritz-certified Arnoldi on the map (right vectors) and on its transpose
    (left), paired by value, each through the dense path's own per-mode body
    (`_floquet_mode`).  `nmodes=None` keeps refusing there: the modal
    spectra need every mode (orbital weight spread over m/n = 0.97 of them),
    and the refusal points to `pnoise`.

    Forced on `_osc_with_ladder(16, 14, 14)` (the instance's limit lowered
    to 8) against the dense modes, gear's pair map (n = 32) and radau's
    (n = 16): multipliers to 7e-15, `p` / `q` to 1e-12.  The certification
    gate is alive: with the Arnoldi budget cut to 12 vectors the four modes
    come back uncertified (Ritz residuals 2e-3 .. 6e-2), warned and flagged
    -- and genuinely wrong, the multipliers 7e-4 .. 0.11 off; at 16 they
    certify (residuals ~1e-20) and are exact."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir, p = _osc_with_ladder(16.0, 14, 14)
    q = PSS(cir, method='radau', reltol=1e-11)
    with quiet():
        q.solve(period=float(p.period), timestep=float(p.period) / 200,
                x0=np.asarray(p._period_state[1], dtype=float)[:cir.n - 1],
                maxiterations=100)
    assert q.converged
    for pp in (p, q):
        with quiet(UsageWarning):
            dense = pp.floquet_modes(nmodes=4)
        pp.FLOQUET_DENSE_LIMIT = 8
        try:
            with _w.catch_warnings(record=True) as rec:
                _w.simplefilter('always')
                ritz = pp.floquet_modes(nmodes=4)
        finally:
            del pp.FLOQUET_DENSE_LIMIT
        assert any('DOMINANT' in str(r.message) for r in rec)
        assert len(ritz) == 4
        for d, r in zip(dense, ritz):
            assert r['certified'] and abs(d['lam'] - r['lam']) < 1e-10, (d['lam'], r['lam'])
            pr, qr = _align_mode(d, r)
            assert np.max(np.abs(pr - d['p'])) < 1e-9 * np.max(np.abs(d['p']))
            assert np.max(np.abs(qr - d['q'])) < 1e-9 * np.max(np.abs(d['q']))
    ## the certification gate, alive: a budget too small to certify
    with quiet(UsageWarning):
        dense = p.floquet_modes(nmodes=4)
    for budget, want in ((12, False), (16, True)):
        p.FLOQUET_DENSE_LIMIT = 8
        p.PPV_RITZ_MAX_BASIS = budget
        try:
            with _w.catch_warnings(record=True) as rec:
                _w.simplefilter('always')
                ritz = p.floquet_modes(nmodes=4)
        finally:
            del p.FLOQUET_DENSE_LIMIT, p.PPV_RITZ_MAX_BASIS
        assert [r['certified'] for r in ritz] == [want] * 4, (budget, ritz)
        assert any('did not certify' in str(r.message) for r in rec) == (not want)
        err = max(abs(d['lam'] - r['lam']) for d, r in zip(dense, ritz))
        assert (err < 1e-10) if want else (err > 1e-4), (budget, err)
    ## the refusals name pnoise for the spectra
    p.FLOQUET_DENSE_LIMIT = 8
    try:
        with pytest.raises(NotImplementedError, match='pnoise'):
            p.floquet_modes()
        with pytest.raises(NotImplementedError, match='pnoise'):
            p.floquet_modes(nmodes=99)
    finally:
        del p.FLOQUET_DENSE_LIMIT


def test_the_ritz_modes_complete_a_conjugate_pair_and_certify_it():
    """`_ritz_modes` never splits a complex-conjugate pair: asked for the
    ONE dominant mode of an operator whose dominant multipliers are a pair,
    it returns both.  Synthetic, because no circuit fixture here has a
    complex dominant multiplier: a 300 x 300 matrix with a known pair at
    0.9 e^{+-0.3 j} over a spectrum inside 0.5."""
    rng = np.random.default_rng(3)
    n = 300
    D = np.diag(0.5 * rng.random(n))
    D[:2, :2] = 0.9 * np.array([[np.cos(0.3), -np.sin(0.3)], [np.sin(0.3), np.cos(0.3)]])
    S = rng.standard_normal((n, n)) / np.sqrt(n) + np.eye(n)
    A = S @ D @ np.linalg.inv(S)
    cir, p = _osc_with_ladder(16.0, 2, 2, npts=40)
    lam, X, res, kk = p._ritz_modes(lambda v: A @ v, n, 1)
    assert len(lam) == 2 and abs(lam[0] - np.conj(lam[1])) < 1e-9, lam
    assert max(res) <= p.PPV_RITZ_RESIDUAL_TOL, res
    ## to the certification level: a Ritz residual of 1e-6 bounds the
    ## eigenvalue error by it (measured 1.1e-7 here)
    for i in range(2):
        assert abs(abs(lam[i]) - 0.9) < 1e-6 and abs(abs(np.angle(lam[i])) - 0.3) < 1e-6
        assert np.linalg.norm(A @ X[:, i] - lam[i] * X[:, i]) < 1e-6


def test_floquet_modes_above_the_dense_limit_on_a_staged_solve_read_the_total_map():
    """On a staged solve the dominant modes are the TOTAL map's (the crossings
    move with the state), on both sides -- as the dense path's
    `total_matrix`.  Forced (the limit lowered to 2) on the comparator
    oscillator against the dense modes."""
    circuit.default_toolkit = circuit.numeric
    osc = _comparator_relaxation_oscillator()
    seed, Tl = _relaxation_oscillator_seed(osc)
    p = PSS(osc, method='radau', reltol=1e-9)
    with quiet(UsageWarning):
        p.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100)
        assert p.converged and p._event_columns is not None
        dense = p.floquet_modes(nmodes=2)
        p.FLOQUET_DENSE_LIMIT = 2
        try:
            ritz = p.floquet_modes(nmodes=2)
        finally:
            del p.FLOQUET_DENSE_LIMIT
    for d, r in zip(dense, ritz):
        assert abs(d['lam'] - r['lam']) < 1e-9, (d['lam'], r['lam'])
        pr, qr = _align_mode(d, r)
        assert np.max(np.abs(qr - d['q'])) < 1e-7 * np.max(np.abs(d['q']))
    assert abs(ritz[0]['lam'] - 1.0) < 1e-6


def test_floquet_modes_run_on_a_monodromy_wider_than_the_dense_limit():
    """Genuinely above `FLOQUET_DENSE_LIMIT`: a van der Pol with a 200-section
    ladder under gear, its pair map 404 wide.  The three dominant modes
    certify, `q^T C p = 1`, and the second multiplier equals `ppv`'s --
    two independent routes (Arnoldi on `M` here, on `I - M` there).
    Measured: 16.5 s to solve, 10.8 s for the modes."""
    circuit.default_toolkit = circuit.numeric
    cir, p = _osc_with_ladder(16.0, 200, 14, npts=60)
    assert p._state_map().width > p.FLOQUET_DENSE_LIMIT
    with quiet(AccuracyWarning, UsageWarning):
        modes = p.floquet_modes(nmodes=3)
        _v, info = p.ppv()
    assert len(modes) == 3 and all(md['certified'] for md in modes)
    assert abs(modes[0]['lam'] - 1.0) < 1e-8
    assert all(md['residual'] < 1e-10 for md in modes)
    lam2 = float(info['second_multiplier'])
    assert abs(abs(modes[1]['lam']) - lam2) < 1e-6, (modes[1]['lam'], lam2)


def test_the_truncated_lam2_is_gated_on_its_own_ritz_residual():
    """Above `FLOQUET_DENSE_LIMIT` a truncated `lam2` must certify itself.

    The dense route closed this below the limit; above it the Arnoldi is all
    there is, and that is exactly where it is least trustworthy — a big circuit
    is the one likely to carry the many slow nodes that break the selection
    (Lai's gated-capacitor DCO: 64 capacitors, >500 equations).

    ⚠ A BIGGER CONSTANT CANNOT BE THE ANSWER, and that is measured, not
    argued: `k = 16` is exact on `_osc_with_ladder(16, 14, 14)` and wrong on
    `(16, 20, 20)` and `(16, 26, 26)`, because the basis has to grow with the
    SLOW-MODE COUNT.  So the basis doubles until the selected pair's own Ritz
    residual `|h_{k+1,k}|·|y_i[last]|/‖y_i‖` certifies it — the same rule at
    every size, and free, since both factors are already in `H`.

    ⚠ IT IS THE EIGENPAIR RESIDUAL, NOT THE SOLVE RESIDUAL.  `_arnoldi_gmres`'s
    own note says a drifted basis "gives multipliers that are wrong in a way
    the residual cannot see" — true of the GMRES residual, false of this one.

    Measured on the truncated path, forced by lowering the dense limit:

        nslow   dense λ₂       gated λ₂       gap ratio   residual   certified
          11    0.995706203    0.995706197      1.000     3.12e-07   yes
          12    0.996324417    0.996324417      1.000     0          yes
          13    0.996818781    0.996818781      1.000     0          yes
          14    0.997220139    0.997220139      1.000     4.0e-76    yes

    and with the basis budget starved so it cannot grow, the recorded failures
    come back **and every one of them is flagged**:

        nslow   gap ratio   residual   certified
          12      -0.031    2.80e-04     no
          13      18.020    3.47e-04     no
          14       0.245    2.08e-03     no

    **Zero false accepts and zero false rejects on this fixture.**  The two
    populations are 3 orders apart here (3.1e-07 against 2.8e-04); a peer's
    independent sweep puts their medians 13 decades apart but ⚠ TOUCHING at
    ~1e-5, which is why `PPV_RITZ_RESIDUAL_TOL` is 1e-6 and not a magic 1e-8.

    ⚠ The residual is a gate on THIS pair, not a truncation bound in general —
    compare `orbital_mode_weights`, whose reconstruction residual saturates at
    whatever the omitted null modes carry.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    fps = {}
    for nslow in (11, 12, 13, 14):
        _cir, pss = _osc_with_ladder(16.0, 14, nslow)
        fp = pss.factored_period()
        ld, _lams = _dense_lam2_of(fp)
        fps[nslow] = (fp, ld, pss)

    def run_all(**settings):
        out = {}
        for nslow in (11, 12, 13, 14):
            fp, ld, pss = fps[nslow]
            ## the settings are this run's alone: each is removed again, so
            ## the next run -- and (5) -- starts from the class defaults
            for name, value in settings.items():
                setattr(pss, name, value)
            try:
                with _w.catch_warnings(record=True) as caught:
                    _w.simplefilter('always')
                    _v, info = pss.ppv()
            finally:
                for name in settings:
                    delattr(pss, name)
            la = float(info['second_multiplier'])
            out[nslow] = dict(
                ratio=(1.0 - la) / (1.0 - ld),
                resid=float(info['second_multiplier_residual']),
                cert=bool(info['second_multiplier_certified']),
                route=info['second_multiplier_route'],
                warned=any('NOT CERTIFIED' in str(c.message) for c in caught))
        return out

    ## force the truncated path on a map the dense route would take
    truncated = dict(FLOQUET_DENSE_LIMIT=4)

    ## (1) WITH ROOM TO GROW: exact everywhere, certified, silent.
    grown = run_all(**truncated)
    for nslow, r in grown.items():
        assert r['route'] == 'arnoldi', \
            'nslow=%d did not reach the truncated path' % nslow
        ## ⚠ A CERTIFIED VALUE IS ACCURATE TO ABOUT ITS RESIDUAL, NOT
        ## TO MACHINE PRECISION, and the two track: nslow=11 certifies at
        ## k=12 with residual 3.1e-07 and lands 1.5e-06 out in the gap.
        ## 12/13/14 come back EXACT because the Krylov space closes on an
        ## invariant subspace there (residual 0), which is a stronger
        ## outcome than the gate promises.
        assert abs(r['ratio'] - 1.0) < 1e-4, \
            'nslow=%d: the gated Arnoldi should track the spectrum and ' \
            'the gap ratio is %.6f. A fixed k=12 gave -0.031 / 18.020 / ' \
            '0.245 at 12/13/14.' % (nslow, r['ratio'])
        assert abs(r['ratio'] - 1.0) < max(1e3 * r['resid'], 1e-9), \
            'nslow=%d: gap error %.2e against a certified residual of ' \
            '%.2e -- the residual is supposed to BOUND the error to ' \
            'within a few orders, and if it stops doing so the gate is ' \
            'certifying something it cannot see' \
            % (nslow, abs(r['ratio'] - 1.0), r['resid'])
        assert r['cert'] and not r['warned'], \
            'nslow=%d: a correct value must certify silently (cert=%s, ' \
            'warned=%s)' % (nslow, r['cert'], r['warned'])

    ## (2) BUDGET STARVED: the recorded failures return, and EVERY one is
    ## flagged.  This is the half that says the gate detects rather than
    ## that the growth happens to help.
    starved = run_all(PPV_RITZ_MAX_BASIS=12, **truncated)
    for nslow, want in ((12, -0.031), (13, 18.020), (14, 0.245)):
        r = starved[nslow]
        assert abs(r['ratio'] - want) < 0.02, \
            'nslow=%d: starved of basis this should reproduce the ' \
            'recorded gap ratio %.3f and gives %.3f' \
            % (nslow, want, r['ratio'])
        assert not r['cert'] and r['warned'], \
            'NO FALSE ACCEPT is the whole claim: nslow=%d is wrong by ' \
            '%.3f in the gap and reported certified=%s / warned=%s' \
            % (nslow, r['ratio'], r['cert'], r['warned'])
    ## and the one that is RIGHT at k=12 must still certify -- a gate that
    ## rejected everything would pass the line above and be useless.
    assert starved[11]['cert'] and abs(starved[11]['ratio'] - 1.0) < 1e-5, \
        'NO FALSE REJECT: nslow=11 is correct at k=12 (residual 3.1e-07) ' \
        'and must still certify; got cert=%s ratio=%.6f' \
        % (starved[11]['cert'], starved[11]['ratio'])

    ## (3) THE RESIDUAL MUST ACTUALLY SEPARATE THEM, or (2) passed for
    ## some other reason.
    ok = starved[11]['resid']
    bad = [starved[k]['resid'] for k in (12, 13, 14)]
    assert ok < PSS.PPV_RITZ_RESIDUAL_TOL < min(bad), \
        'the residual no longer brackets the tolerance: right %.2e, ' \
        'wrong %r, tol %.0e' % (ok, bad, PSS.PPV_RITZ_RESIDUAL_TOL)
    assert min(bad) / ok > 100.0, \
        'right and wrong are only %.1fx apart in residual (%.2e vs %.2e); ' \
        'with that little margin the tolerance is a tuned constant rather ' \
        'than a separation' % (min(bad) / ok, ok, min(bad))

    ## (4) ⚠ NEUTER THE TOLERANCE: with it wide open the wrong values must
    ## certify, which is what proves the tolerance is load-bearing and not
    ## decoration.
    blind = run_all(PPV_RITZ_RESIDUAL_TOL=1.0, PPV_RITZ_MAX_BASIS=12,
                    **truncated)
    assert all(blind[k]['cert'] for k in (12, 13, 14)), \
        'with the tolerance opened to 1.0 the wrong values should sail ' \
        'through; if they do not, something other than the tolerance is ' \
        'gating and this test is not measuring what it says'
    assert abs(blind[14]['ratio'] - 0.245) < 0.02, \
        'and they should be the SAME wrong values (%.3f)' \
        % blind[14]['ratio']

    ## (5) AND THE DEFAULT PATH IS UNTOUCHED: n = 32 is inside the real limit,
    ## so this fixture still takes the spectrum and certifies trivially.
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        _v, info = fps[14][2].ppv()
    assert info['second_multiplier_route'] == 'dense' \
        and info['second_multiplier_certified'] \
        and info['second_multiplier_residual'] == 0.0, \
        'the dense route must report itself exact and certified, got %r' \
        % {k: info[k] for k in ('second_multiplier_route',
                                'second_multiplier_certified',
                                'second_multiplier_residual')}
    assert not any('NOT CERTIFIED' in str(c.message) for c in caught), \
        'the dense route must not warn'


def test_the_twin_serves_the_whole_frequency_aware_ppv():
    """On trap and euler the monodromy TWIN serves `ppv()`, `floquet_modes`
    and `factored_period()` -- and so, whole, `frequency_aware_ppv`: a trap
    solve's frequency-aware PPV is its twin's, bit for bit.

    ⚠ Until 2026-09-27 it borrowed only the twin's period map and ran the
    propagation on the trap solve itself: exact on an ODE, but with an
    ALGEBRAIC node (here the realisation `white_ref`, a source behind a
    multiplier) the samples read 1.6e-3 off the twin's own, and `c(f)`
    1.3e-4 -- varying with the offset, so no carrier-frequency factor.
    """
    cir, pss, pac, ov = _orbit_modulated_vdp('white_ref', method='trap')
    tw = pss.monodromy_twin()
    assert tw is not pss
    f = 1e-2 / float(pss.period)
    with quiet():
        _va, ia = pss.frequency_aware_ppv(f)
        _vb, ib = tw.frequency_aware_ppv(f)
    assert np.array_equal(np.asarray(ia['samples_eq']), np.asarray(ib['samples_eq']))


def test_floquet_modes_runs_under_the_stage_methods_and_conserves_qCp():
    """`floquet_modes` CRASHED under trbdf2 for as long as the DIRK transposed
    matvec existed -- "can't multiply sequence by non-int of type 'complex'".

    The complex-vector `collect` path recombined real and imaginary parts with
    a flat `a + 1j*b` over a NESTED list of per-stage solves (with `None` for
    the explicit first stage).  The same line sat at FOUR sites, one per
    transposed-matvec variant.  Nothing exercised it: every Floquet test in
    this file used gear.  Fixed 2026-09-07 with one recursive `_cx_collect`.

    This asserts the two stage methods produce modes at all, and that the
    C-weighted invariant holds at the same index -- the property the `C^-T`
    fix of the same day restored.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    for method in ('trbdf2', 'radau'):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
        pss = PSS(cir, method=method, reltol=1e-12)
        with quiet():
            pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                      maxiterations=300)
            modes = pss.floquet_modes()      # used to raise under trbdf2
        assert len(modes) == 2, '%s: expected 2 modes, got %d' % (
            method, len(modes))
        x0r = np.delete(np.asarray(pss.waveform[1], dtype=float)[:, 0],
                        pss.irefnode)
        Cm = np.asarray(pss._C_at(x0r), dtype=float)
        for k, md in enumerate(modes):
            n_ = min(md['p'].shape[1], md['q'].shape[1])
            cyc = [abs(complex(np.vdot(md['q'][:, j], Cm @ md['p'][:, j])) - 1.0)
                   for j in range(0, n_ - 1, max(1, n_ // 8))]
            assert max(cyc) < 5e-2, \
                '%s mode %d: q^T C p drifts by %.3e around the cycle' % (
                    method, k, max(cyc))


def test_injection_locking_range_on_the_hostile_fixture_is_set_by_its_ppv_fundamental():
    """A6, second gate: the hostile fixture locks over TWICE the averaging
    range, and the reason is measured by an independent route.

    Same instrument as the control gate.  On the hostile fixture (van der Pol
    + 0.30 u^2, C = 4, L = 1/4) the locked branch is tracked to 2.00x the
    averaging prediction I_inj/(2CA) and lost at 2.10x (measured 2026-09-07),
    where the control and a SYMMETRIC C = 4 fixture both give exactly 1.00x.
    So the factor is the ASYMMETRY, not a C-dependence in the prediction --
    the symmetric C = 4 discriminator is what ruled that out, and it was run
    because a clean factor of two on a fixture the control cannot distinguish
    is a warning, not a result.

    ⚠⚠ THE MECHANISM, BY AN INDEPENDENT ROUTE.  The lock range is the
    injection times the PPV's FUNDAMENTAL per unit current, |Gamma_1|/C:

        control C=1        0.2500 / 1  =  0.0625 ... = 1/(2CA)   exactly
        symmetric C=4      0.2500 / 4  =  0.0625     = 1/(2CA)   exactly
        HOSTILE C=4        0.4588 / 4  =  0.1147     = 1.83x

    Four digits on both symmetric fixtures with no fit -- the averaging
    result IS the PPV fundamental, tying A6 to A2.  On the hostile fixture the
    PPV predicts 1.83x against the sweep's 2.0-2.1x: a ~10% gap, and the
    hostile PPV's SECOND harmonic is 10% of its fundamental (|G2|/|G1| = 0.10,
    against 1e-4 on both symmetric orbits).  A first-order phase prediction on
    an orbit with that much second-harmonic sensitivity should miss by about
    that much.  Attributed, not proven; the bound below admits it.
    """
    last, pred_hz, f0, trace = _injection_lock_edge(
        _vdp_injected(cval=4.0, lval=0.25, a=0.30),
        steps=np.arange(0.0, 2.61, 0.2))
    assert last is not None, 'no locked branch found: %r' % (trace[:3],)
    assert 1.8 <= last <= 2.2, \
        'hostile locked branch tracked to %.2fx the averaging prediction; ' \
        'measured 2.00x (lost at 2.10x). Trace: %r' % (
            last, [(r, round(am, 3), round(lm, 5), lk) for r, am, lm, lk in trace])

    ## the PPV fundamental, per unit current, on the same three fixtures

    def gamma1(cval, lval, a):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=cval)
        cir['L'] = L('v', gnd, L=lval)
        mu = 1.0 / (2.0 * np.pi * 8.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + a * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        T0 = 2.0 * np.pi * np.sqrt(cval * lval)
        p = PSS(cir, method='gear', reltol=1e-11)
        with quiet(AccuracyWarning):
            p.solve(period=T0, timestep=T0 / 400, x0=np.array([2.0, 0.0]),
                    maxiterations=400)
            v0, info = p.ppv()
        iv = cir.get_node_index('v')
        ivr = iv if iv < p.irefnode else iv - 1
        g = np.array([np.asarray(v0, float)[ivr]]
                     + [np.asarray(s_, float)[ivr] for s_ in info['samples']])[:-1]
        G = np.fft.rfft(g) / len(g)
        W = np.asarray(p.waveform[1], float)[iv]
        A = (W.max() - W.min()) / 2
        return abs(G[1]) / cval, abs(G[2]) / abs(G[1]), 1.0 / (2.0 * cval * A)

    g_c, h_c, pred_c = gamma1(1.0, 1.0, 0.0)
    g_s, h_s, pred_s = gamma1(4.0, 0.25, 0.0)
    g_h, h_h, pred_h = gamma1(4.0, 0.25, 0.30)
    ## 1. on BOTH symmetric fixtures the PPV fundamental IS 1/(2CA)
    for lbl, g, pr in (('control', g_c, pred_c), ('symmetric C=4', g_s, pred_s)):
        assert abs(g / pr - 1.0) < 2e-3, \
            '%s: |Gamma_1|/C = %.6f against 1/(2CA) = %.6f' % (lbl, g, pr)
    ## 2. on the hostile fixture it is ~1.8x, and that is the lock-range factor
    ratio_ppv = g_h / pred_h
    assert 1.7 <= ratio_ppv <= 2.0, \
        'hostile PPV fundamental is %.3fx 1/(2CA); measured 1.83x' % ratio_ppv
    ## 3. and the asymmetry is visible where it should be: a second-harmonic
    ##    PPV component the symmetric orbits do not have
    assert h_h > 0.05 and h_c < 1e-3 and h_s < 1e-3, \
        '|Gamma_2|/|Gamma_1|: control %.1e, symmetric C=4 %.1e, hostile %.3f' % (
            h_c, h_s, h_h)


def test_the_frequency_aware_equation_rows_are_batched_bit_for_bit():
    """`frequency_aware_ppv`'s equation-row samples (2026-09-27): the
    per-sample blocks of `C(x_j)` / `G(x_j)` are read once per solved
    orbit and the solves run stacked (`_equation_row_blocks`,
    `_equation_row_batch`), instead of re-evaluating the circuit at every
    sample of every solve (36 % of a solve under the profiler, ~10 % in
    wall time; with the lazy DC line, the all-orders lineshape 14.4 ->
    9.5 s).  BIT-identical to the per-sample path at four offsets.  The
    circuit makes every block move along the orbit: a tanh charge at `v`
    (`C` varies by 0.27) and a cubic conductance at the ALGEBRAIC node `x`
    (its fill's `G` block varies by 0.20).  On the plain slow-node LC they
    are constant, and a batch reading the wrong sample's blocks passed.
    The re-solve test above covers the cache following a new orbit."""
    from pycircuit.circuit.shooting import PSS as _PSS
    circuit.default_toolkit = circuit.numeric
    T0 = 6.6634
    c = SubCircuit()
    c.add_node('v'); c.add_node('w'); c.add_node('x')
    c['C'] = C('v', gnd, c=1.0)
    c['Q'] = BSource('v', gnd, gnd, 'v', q_func=lambda u: 0.3 * np.tanh(u))
    c['L'] = L('v', 'x', L=1.0)
    c['Rl'] = R('x', gnd, r=0.2)
    c['Gx'] = BSource('x', gnd, gnd, 'x', i_func=lambda u: 0.5 * u ** 3)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: (u - u ** 3 / 3.0) + 0.25 * (u ** 2 - 2.0))
    c['Rs'] = R('v', 'w', r=1e2, noisy=False)
    c['Cs'] = C('w', gnd, c=100.0 * T0 / 1e2)
    c['nw'] = IS('w', gnd, i=0.0, noisePSD=1e-6)
    pss = PSS(c, method='gear', reltol=1e-11)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=T0, timestep=T0 / 240, x0=x0, maxiterations=200)
    assert pss.converged
    f0 = 1.0 / float(pss.period)
    with quiet(AccuracyWarning):
        for nu in (1e-5, 1e-3, 3e-2, 0.4):
            a = pss.frequency_aware_ppv(nu * f0)[1]['samples_eq']
            blocks = pss._eq_row_cache[2]
            assert blocks is not None, 'the batch did not run'
            assert blocks['A'] is not None and blocks['rows'], 'no algebraic fill'
            for k in ('C', 'A'):
                assert np.ptp(blocks[k], axis=0).max() > 0.1, \
                    '%s does not move along the orbit' % k
            orig = _PSS._equation_row_blocks
            unbatched = []
            _PSS._equation_row_blocks = lambda self, *args: unbatched.append(1)
            try:
                b = pss.frequency_aware_ppv(nu * f0)[1]['samples_eq']
            finally:
                _PSS._equation_row_blocks = orig
            ## ⚠ vacuous if the patch is never looked up
            assert unbatched, 'the patched _equation_row_blocks was never consulted'
            assert np.array_equal(a, b), (nu, float(np.max(np.abs(a - b))))


def test_the_continuous_adjoints_arnoldi_path_equals_its_dense_path():
    """Item 3 of gear-first-class-on-non-uniform-grids: above a small `m` the
    continuous adjoint no longer forms the dense 2m x 2m backward map but runs
    Arnoldi on it, Ritz-certified.  Forced onto the hostile fixture (m = 2) by
    lowering the switch, the Arnoldi q must equal the dense q -- MEASURED cos
    1.00000000 and identical invariant spreads at a basis of 3 / 16 / 24 for
    2m = 4 / 32 / 124 on the hostile and ladder fixtures.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)

    def fracs(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=4.0)
    cir['L'] = L('v', gnd, L=0.25)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300, break_events=False, grid=fracs(400))
    assert pss.converged
    pss.CONTINUOUS_ADJOINT_DENSE_M = 8
    dense = pss.floquet_modes()
    pss.CONTINUOUS_ADJOINT_DENSE_M = 0          # force Arnoldi at m = 2
    arn = pss.floquet_modes()
    assert len(dense) == len(arn) >= 2
    for a, b in zip(dense, arn):
        qa, qb = np.asarray(a['q']), np.asarray(b['q'])
        for j in (0, 100, 200, 399):
            cos = abs(np.vdot(qa[:, j], qb[:, j])) / (np.linalg.norm(qa[:, j]) * np.linalg.norm(qb[:, j]))
            assert cos > 1.0 - 1e-8, (j, cos)
        assert np.linalg.norm(qa - qb) / np.linalg.norm(qa) < 1e-6


def test_gear_runs_its_ppv_on_the_lte_grid_gear_produced_and_is_told_when_its_unit_multiplier_left_the_circle():
    """Gear on an ADAPTIVE grid (2026-09-21, Andreas: "Do as you suggest").
    `lte_grid`'s adaptive run is `Gear2Integrator()` by default, so the grid
    is gear-shaped -- and gear could not run its PPV on it.

    Relaxation van der Pol, mu = 10, `lte_grid` at reltol 1e-5, ~195 steps,
    span 176:1.  Two window defects: the run's FINAL step is truncated to
    land on `tend` (a tenth of its predecessor), so a periodic grid carried
    a 10.5x growth across the seam -- beyond the two-step zero-stability
    bound -- and the integrator dropped two steps to Euler, which gear's
    transposed replay refuses; rotating the seam only moved that step into
    the interior.  Fixed in `lte_grid`: the window ends at the last natural
    step and, when the seam is outside the bound, the cut is rotated by the
    least to a sane one (an unconditional move to the flattest seam put the
    seed at a phase from which the mu = 4 test's Newton, 8 % off in period,
    ran away).  Two more things stood between gear and its own grid:
    `_period_grid`'s first-step subdivision put ONE step of `fr.min()` in
    front of a coarse first step (a 138x growth) -- now a doubling ramp, every
    ratio 2 -- and Transient's `h_curr/h_last < 0.1` drop to Euler (a
    stalled-estimate heuristic for adaptive runs) fired on the ramp's first
    step at the seam -- `Gear2Integrator.shrink_guard`, off in PSS-driven
    transients.  Measured after: seam inside the bound, max growth 2.000, no
    two-alpha steps, and gear's PPV, modes and c run.

    ⚠ THE PERIOD PASSED MUST BE CLOSE.  With it 16 % off, gear's AND
    trbdf2's per-step Newton fail on the coarse steps that land on the
    edges (the grid is fractions of the period); seeded at the reference
    both converge.  Measured, 195 points, against radau uniform N = 1480
    (self-checked against 740 to 5e-3 ppm)::

        method / grid          T err (ppm)   c rel     |rho - 1|
        gear   LTE             -519          +0.52     1.0e-1
        gear   LTE split 2x    -133          +0.15     2.7e-2
        gear   uniform, 185    -7214         +0.10     --
        trbdf2 LTE             -66           +0.16     1.4e-2
        trbdf2 LTE split 2x    -17           +0.040    3.4e-3

    Gear's period is second order on the adaptive grid and 14x better than
    a uniform grid of the same count; its c is 52 % high because the
    discrete period map's unit multiplier sits 0.10 off the circle, which
    `ppv` -- solved at exactly 1 -- could not see: it now WARNS from
    `spectral_radius` (`PPV_UNIT_MULTIPLIER_WARN`), c tracking the
    departure at 5-12x.  The spline quadrature is not it (trapezoid +0.528).
    And trbdf2 uses gear's own grid 8x better than gear for the period, 3.5x
    for c: on a relaxation oscillator's adaptive grid, that is the method.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    MU = 10.0
    T_REF = 19.098600502

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: MU * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    cir = vdp()
    p = PSS(cir, method='gear')
    xfull = np.zeros(cir.n)
    xfull[[str(n_) for n_ in cir.nodes].index('v')] = 2.0
    with quiet():
        fr, seed = p.lte_grid(T_REF, x0=xfull, reltol=1e-5)
    fr = np.asarray(fr, float)
    N = len(fr)
    r = fr[1:] / fr[:-1]
    ## ⚠ 2026-09-28 `lte_grid` designs on its own `relref`, default
    ## 'pointlocal': 273 points (205 on 'sigglobal', which gear's grid had
    ## used since D3), and better -- gear -6.9 ppm (2x +3.3), rho - 1 0.041,
    ## trbdf2 +4.1 ppm, against +17.9 (+12.5), 0.061, +10.4 on 'sigglobal'
    assert 150 < N < 340, N
    assert 1.0 / 2.414 < fr[0] / fr[-1] < 2.414, fr[0] / fr[-1]     # was 10.53
    assert r.max() < 2.5, r.max()                                     # the controller's clamp

    def run(method, grid, n):
        c_ = vdp()
        q = PSS(c_, method=method, reltol=1e-9)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=T_REF, timestep=T_REF / n, x0=seed, maxiterations=80,
                    break_events=False, grid=grid)
            assert q.converged
            fp = q.factored_period()
            if fp.kind == 'solved_history':
                assert all(len(st[2]) == 3 for st in fp.steps), 'an order-dropped step'
            q.ppv()
        warned = [str(w_.message) for w_ in rec if 'unit multiplier sits' in str(w_.message)]
        return 1e6 * (float(q.period) - T_REF) / T_REF, float(q.spectral_radius), warned
    e1, rho1, w1 = run('gear', fr, N)
    e2, rho2, w2 = run('gear', np.repeat(fr / 2.0, 2), 2 * N)
    ## on the FOLDED grid (every settled period, per-phase minimum, seam at
    ## the finest phase; 2026-09-21 "Do 1"): gear +20 / +12 / +4.4 ppm under
    ## 1x / 2x / 4x -- against -1461 / -22 on two raw single windows, a
    ## uniform grid of the same count at -6541 -- and its unit multiplier
    ## nearer the circle (c +22 %, was +52 %); still told
    ## gear's period on a ~200-point folded grid is an O(300 ppm) constant
    ## whose sign depends on the remesh phase (+285 with the seam on the
    ## slow branch, +20 / -28 elsewhere), second order under splitting
    assert abs(e1) < 500 and abs(e2) < 150 and abs(e2) < abs(e1), (e1, e2)
    ## and it is told exactly when its unit multiplier has left the circle
    for rho, w in ((rho1, w1), (rho2, w2)):
        assert bool(w) == (rho - 1.0 > PSS.PPV_UNIT_MULTIPLIER_WARN), (rho, w)
    e3, rho3, w3 = run('trbdf2', fr, N)
    assert abs(e3) < 100, (e3, e1)                                  # +7 .. +57 by remesh phase
    assert bool(w3) == (rho3 - 1.0 > PSS.PPV_UNIT_MULTIPLIER_WARN), (rho3, w3)


def test_a_staged_glm_oscillators_ppv_and_floquet_modes_carry_the_moving_events():
    """A Nordsieck GLM oscillator with state events: `ppv()` and
    `floquet_modes` RAISED ("costate injections are not built on a
    Nordsieck GLM (it has no state-event stage)", false by then), because
    the GLM's own map on the state (`_GLMStateMap`) refused the event rows'
    costate injections.  Two rounds had each been sound alone: the GLM's
    event stage, and the GLM reading its own map for `ppv`.  Now the
    injection at node n seeds step n-1's last stage (which IS `x_n`), node
    0's the result (`_GLMPeriod.x_matvec_transposed`).  Measured against
    the EXACT saltation PPV (`_exact_relaxation_oscillator_ppv`) at 200
    points: glm3 3.9e-4, glm2 6.6e-4, radau 4.5e-4 -- the comparison's own
    floor (see the radau test below); the second multiplier 0.02 for all.
    """
    circuit.default_toolkit = circuit.numeric
    T_ex, exact = _exact_relaxation_oscillator_ppv()
    for method in ('glm3', 'glm2'):
        cir = _comparator_relaxation_oscillator()
        names = [str(n_) for n_ in cir.nodes]
        seed, Tl = _relaxation_oscillator_seed(cir)
        q = PSS(_comparator_relaxation_oscillator(), method=method, reltol=1e-9)
        with quiet(UsageWarning):
            q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100,
                    state_events=True)
            v, info = q.ppv()
            fm = q.floquet_modes(nmodes=2)
        assert q.converged and len(q._state_event_fracs) == 4, method
        red = [n_ for i, n_ in enumerate(names) if i != q.irefnode]
        idx = [red.index(nm) for nm in ('c', 'fb0', 'fb1')]
        X = np.asarray(q.waveform[1], dtype=float)
        ts = np.asarray(q.waveform[0], dtype=float)
        S = np.asarray(info['samples'])
        nodes = list(q._event_columns['nodes'])
        before = [j for j in (5, 15, 25) if j < nodes[0]]
        after = [int(np.searchsorted(ts, f * float(q.period))) for f in (0.5, 0.8)]
        worst = 0.0
        for j in before + [nodes[1] + 3] + after:
            ex = exact(X[[names.index(nm) for nm in ('c', 'fb0', 'fb1')], j])
            worst = max(worst, float(np.max(np.abs(S[j, idx] - ex))
                                     / np.max(np.abs(ex))))
        assert worst < 2e-3, (method, worst)
        lam = sorted((abs(x['lam']) for x in fm), reverse=True)
        assert abs(lam[0] - 1.0) < 1e-6 and abs(lam[1] - 0.02) < 1e-3, (method, lam)
        ## the injection is EXACT, by duality with the forward pass: the
        ## transposed map with `inject` is the gradient of `v . x_N + sum_n
        ## inject[n] . x_n` (on this fixture the events' weight in the PPV
        ## is small, so accuracy against the exact PPV cannot see it)
        sm = q.factored_period().state_map()
        m = sm.width
        rng = np.random.default_rng(3)
        vv, u = rng.standard_normal(m), rng.standard_normal(m)
        inj = rng.standard_normal((len(sm.steps), m))
        grad = np.asarray(sm.matvec_transposed(vv, inject=inj))
        _xN, fwd = sm.forward_states(u)
        lin = float(vv @ fwd[-1]) + float(inj[0] @ u) + sum(
            float(inj[n] @ fwd[n - 1]) for n in range(1, len(sm.steps)))
        assert abs(float(grad @ u) - lin) < 1e-10 * abs(lin), (method, float(grad @ u), lin)


def test_the_ppv_on_a_staged_autonomous_solve_is_bordered_and_matches_the_exact_saltation_ppv():
    """Events phase B (2026-09-22): `ppv()` on a staged solve takes the
    null vector of the TOTAL monodromy and carries the crossings' motion
    into its samples as costate injections at the event nodes.  The
    reference is EXACT -- `_exact_relaxation_oscillator_ppv`, the
    piecewise-linear flow with saltation matrices -- and independent of
    every integrator.  Measured with VSwitch's compact transition
    (2026-09-22): bordered and fixed-grid samples both within 4.7e-4 of
    it at 200 points and 2.9e-4 at 400 (the nearest-state comparison's
    floor; INSIDE a 0.04 ns window the discrete sample interpolates the
    jump and the comparison is meaningless).  With the tanh the fixed-grid
    PPV of the same staged solve was 140-230 % off on fb1 before the
    first crossing (`M^T v - v` = 1.5, |lambda - 1| = 0.99): that was the
    transition's tails outside the window, which the bordering corrected
    and the compact transition removed.  ⚠ The transient FD that was to be the
    instrument read `c` at +4.4e-8 for the exact -7.4e-9 s/V -- its phase
    was read at the `c` waveform's crossing of `(max + min) / 2` of the
    PERTURBED record, a level the perturbation moved; read at the
    comparator's own crossing it agrees with the exact value to 0.4 %.
    Also pins the autonomous stage listing its crossings in
    `event_times`."""
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    names = [str(n_) for n_ in cir.nodes]
    seed, Tl = _relaxation_oscillator_seed(cir)
    c2 = _comparator_relaxation_oscillator()
    q = PSS(c2, method='radau', reltol=1e-9)
    with quiet():
        q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100, state_events=True)
    assert q.converged
    T_ex, exact = _exact_relaxation_oscillator_ppv()
    assert abs(q.period / T_ex - 1.0) < 5e-4, (q.period, T_ex)
    ## the autonomous stage lists its four landed crossings
    ev = [float(e) for e in q.event_times]
    assert len(ev) == 4 and all(0.0 < e < 1.0 for e in ev), ev
    red = [n_ for i, n_ in enumerate(names) if i != q.irefnode]
    idx = [red.index(nm) for nm in ('c', 'fb0', 'fb1')]
    X = np.asarray(q.waveform[1], dtype=float)
    ts = np.asarray(q.waveform[0], dtype=float)
    with quiet():
        v, info = q.ppv()
    S = np.asarray(info['samples'])
    nodes = list(q._event_columns['nodes'])
    before = [j for j in (5, 15, 25) if j < nodes[0]]
    after = [int(np.searchsorted(ts, f * float(q.period))) for f in (0.5, 0.8)]
    worst = 0.0
    for j in before + [nodes[1] + 3] + after:
        xs = X[[names.index(nm) for nm in ('c', 'fb0', 'fb1')], j]
        ex = exact(xs)
        got = S[j, idx]
        worst = max(worst, float(np.max(np.abs(got - ex)) / np.max(np.abs(ex))))
    assert worst < 2e-2, worst
    ## the fixed-grid PPV of the same solve is O(1) wrong before the crossings
    cols = q._event_columns
    q._event_columns = None
    try:
        with quiet():
            _vu, info_u = q.ppv()
    finally:
        q._event_columns = cols
    Su = np.asarray(info_u['samples'])
    ## ⚠ RE-PINNED 2026-09-22 for VSwitch's COMPACT transition: with the whole
    ## transition inside the landed window the fixed-grid map of the staged
    ## solve carries the switching itself -- its unit multiplier sits at
    ## 1e-7 / 7e-10 / 7e-10 (100 / 200 / 400 points, the total's at 1e-12 ..
    ## 1e-15) and its PPV matches the exact one as well as the bordered does
    ## (4.7e-4 both at 200 points, the nearest-state comparison's floor).
    ## The 140 % it was off with the tanh was the tails outside the window.
    ## Pinned: the fixed-grid PPV within 2e-2 too, the bordered within 1e-6
    ## of it at the period node (the two null vectors agree to the solve).
    worst_u = 0.0
    for j in before + after:
        xs = X[[names.index(nm) for nm in ('c', 'fb0', 'fb1')], j]
        ex = exact(xs)
        worst_u = max(worst_u, float(np.max(np.abs(Su[j, idx] - ex)) / np.max(np.abs(ex))))
    assert worst_u < 2e-2, worst_u


def test_the_monodromy_twin_is_capped_and_warns_when_it_hits_the_cap():
    """`PSS.TWIN_MAXITER` (Andreas, 2026-09-23): a monodromy twin is a polish
    from this run's converged state, not a cold solve, and it used to inherit
    the caller's `maxiterations` -- on B16's van der Pol a poor (Euler) seed
    then spent 790 s running two twins to their 300-iteration budgets to
    refuse.  The twin is now capped at `min(maxiterations, TWIN_MAXITER)`, 40
    by default, and a twin that hits the cap without converging WARNS before
    it refuses.  Pinned both ways: from a GOOD seed (trapezoidal, the solve
    allowing 300) the capped twin converges -- in 4-13 iterations, measured --
    with no cap warning; from the poor Euler seed, with the cap patched to 6
    under a solve allowing 20, the twin warns that it hit the cap and then
    refuses."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)

    def build():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u ** 2)
        return cir

    def solve(method, maxiterations):
        pss = PSS(build(), method=method, reltol=1e-12)
        with quiet(AccuracyWarning):
            pss.solve(period=6.731, timestep=6.731 / 400,
                      x0=np.array([2.0, 0.0]), maxiterations=maxiterations)
        assert pss.converged
        return pss

    assert PSS.TWIN_MAXITER == 40
    good = solve('trap', 300)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        _v, info = good.ppv()
    assert info['monodromy_method'] == 'radau'
    assert not [r for r in rec if 'capped iteration budget' in str(r.message)]
    poor = solve('euler', 20)
    poor.TWIN_MAXITER = 6
    with pytest.warns(RuntimeWarning, match='within its capped iteration budget'):
        with pytest.raises(RuntimeError, match='did not converge|spurious|too poor'):
            poor.ppv()


def test_radau_and_esdirk43_serve_as_the_monodromy_twin():
    """A trap or euler oscillator reads its monodromy off a TWIN; the twin
    could be trbdf2 or gear only.  Any method whose own map serves --
    `carries_own_monodromy` -- may be one now.  Measured on van der Pol, a
    trap run at 200 points against radau at 800: lambda2 / c / the
    oscillator covariance's d off by 1e-3 / 1e-4 / 2e-3 under a gear twin,
    1e-4 / 6e-6 / 8e-5 under trbdf2, 1e-10 / 2e-11 / 6e-13 under radau at
    about trbdf2's cost.  A Nordsieck GLM (its map on the Nordsieck state)
    is refused.
    """
    circuit.default_toolkit = circuit.numeric
    lam2_ref = 0.20038546770896304          # radau, 800 points
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    p = PSS(_vdp_asym(), method='trap', reltol=1e-10)
    p.monodromy = 'radau'
    with quiet(AccuracyWarning):
        p.solve(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]),
                maxiterations=60)
        _v, info = p.ppv()
    assert info['monodromy_method'] == 'radau'
    assert abs(float(info['second_multiplier']) / lam2_ref - 1.0) < 1e-8
    p.monodromy = 'glm2'
    with pytest.raises(ValueError, match="'radau'"):
        p.monodromy_twin()


def test_the_ppv_reads_a_stateful_diode_at_each_orbit_point():
    """The PPV's own device reads (the algebraic rows' `G`, `q`, `C(0)`)
    take the devices AT the point, as the map's reads do (`devices_at`).
    A van der Pol tank clamped by a diode on an ALGEBRAIC node: while the
    diode is reverse-biased, current into that node can only flow back
    through the resistor into the tank, so its equation-row PPV equals the
    tank's (`v_b = v_v / (1 + R g_d)`, `g_d ~ 0`).  ⚠ Measured before
    (2026-09-28): `Diode.G` linearised at the voltage the last Newton left
    (conducting, at the period's end) on every sample -- the node's entry
    98 % off there, and the diffusion constant 3.4x too large."""

    from pycircuit.circuit.elements import BSource, Diode
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir.add_node('b')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: u - u ** 3 / 3.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['R1'] = R('v', 'b', r=10.0)
    cir['D'] = Diode('b', gnd)
    pss = PSS(cir, method='radau', reltol=1e-10)
    names = [str(n) for n in cir.nodes if str(n) != 'gnd!']
    x0 = np.zeros(cir.n - 1)
    x0[names.index('v')] = 2.0
    with quiet():
        pss.solve(period=2 * np.pi, timestep=2 * np.pi / 400, x0=x0,
                  maxiterations=60)
        assert pss.converged
        _v, info = pss.ppv()
    S = np.asarray(info['samples_eq'], dtype=float)
    vb = np.asarray(pss.waveform[1], dtype=float)[
        [str(n) for n in cir.nodes].index('b')][:len(S)]
    off = vb < 0.3                                   # the diode reverse-biased
    assert off.sum() > 100
    iv, ib = names.index('v'), names.index('b')
    err = (np.max(np.abs(S[off, ib] - S[off, iv]))
           / np.max(np.abs(S[off, iv])))
    assert err < 1e-4, err
