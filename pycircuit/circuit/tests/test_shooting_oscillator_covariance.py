"""Shooting tests: shooting oscillator covariance.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.shooting._noise_components import psd_sqrt
from pycircuit.circuit.tests._shooting_fixtures import (_Flicker,
    _comparator_relaxation_oscillator,
    _exact_relaxation_oscillator_model,
    _loss_osc,
    _orbit_modulated_vdp,
    _rc_noisy,
    _relaxation_oscillator_seed,
    _sampler_fixture,
    _solve_vdp_noise,
    _sw,
    _vdp_asym)


def _osc_cov(npts=240):
    _cir, pss, pac = _solve_vdp_noise(npts=npts)
    ## (pair=True: these tests walk and project in the integrator's pair
    ## space; the default is the m x m node covariance since 2026-09-29)
    K_orb, info = pac.oscillator_covariance(pss, samples=True, pair=True)
    d = info['d']
    return pss, pac, K_orb, d, info


def test_the_oscillator_covariance_predicts_the_walk_forty_periods_out():
    """⚠ THE ASSERTION IS A PREDICTION, NOT A PROPERTY OF THE SOLVE.

    `covariance` refuses an oscillator because `I − M⊗M` is singular there
    and no periodic covariance exists. The claim this makes instead is a
    split — a bounded part plus a random walk along the orbit:

        K(t₀ + nT) = K_orb + n·d·u·uᵀ

    which is falsifiable in a way "the residual is small" is not: run the
    real Lyapunov recursion forward for forty periods (9,600 steps, from
    `K = 0`, touching nothing the bordered solve produced) and compare.
    A bordered system can always be solved; only this says the answer
    means anything.

    ⚠ THE FIRST PERIOD IS THE LOOSEST, AND THAT IS PHYSICS. Starting from
    `K = 0` rather than `K_orb` leaves a transient in the bounded part; it
    decays with `|λ₂| = 8.5e-4` and is gone by period two. Measured
    2.6e-07, 1.1e-10, 2.2e-10, 5.4e-10, 1.2e-09, 2.5e-09 at periods
    1/2/5/10/20/40 — the walk dominates 99.9% of the trace by then, so the
    late numbers test the growth term and the early one tests the bound.
    """
    pss, pac, K_orb, d, info = _osc_cov()
    u = info['tangent_pair']
    As, Qs, _K1, M, _m, n = pac._lyapunov_pieces(pss, 'test')
    ## the split's own periodicity statement: P is periodic UP TO the
    ## growth, which appears exactly once per period -- not periodic
    P = info['samples']
    grow = d * np.outer(u, u)
    assert np.max(np.abs(P[-1] - K_orb - grow)) < 1e-10 * np.max(np.abs(grow))

    K = np.zeros((n, n))
    seen = {}
    for p in range(1, 41):
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        if p in (1, 2, 40):
            pred = K_orb + p * grow
            seen[p] = np.max(np.abs(K - pred)) / np.max(np.abs(pred))
    assert seen[1] < 1e-5, \
        'one period off by %.3e; the bounded part is wrong' % seen[1]
    assert seen[40] < 1e-7, \
        'forty periods off by %.3e; the growth RATE is wrong -- that ' \
        'error accumulates where the first-period one does not' % seen[40]
    assert seen[2] < seen[1], \
        'the initial transient does not decay, so the split is not a ' \
        'decomposition into a settling part plus a walk'


def test_the_growth_rate_is_the_diffusion_constant_by_another_route():
    """⚠ TWO SEPARATELY ANCHORED QUANTITIES, CLOSED INTO A LOOP.

    A phase deviation `α` displaces the state by `α·u`, so the growing
    covariance is `Var(α)·u uᵀ = c·t·u uᵀ` and therefore `d = c·T`. The two
    sides share the `CY/2` convention and nothing else: `c` is a quadratic
    form in the ADJOINT-replayed PPV, `d` comes from a FORWARD Lyapunov
    recursion closed by a bordered Kronecker solve.

    ⚠ AND THEIR ANCHORS ARE INDEPENDENT TOO, which is the property that
    was missing when a 2× error survived a 0.9965 agreement. The
    injection behind `d` is pinned by `kT/C`; `c` is pinned by a nonlinear
    Monte Carlo reading phase from zero crossings. Neither measurement can
    influence the other's answer.

    ⚠ ASSERTED AS A CONVERGENCE, NOT A TOLERANCE. Both carry an O(h)
    piecewise-constant approximation to white noise, so at any single grid
    they differ by a real amount; what must hold is that the difference
    HALVES per doubling. Measured 1.87% → 1.03% → 0.54%. A shared error
    would cancel and give a flat ratio of 1; a wrong one would not shrink.
    """
    errs = []
    for npts in (120, 240, 480):
        _cir, pss, pac = _solve_vdp_noise(npts=npts)
        _K, info = pac.oscillator_covariance(pss)
        c = pac.diffusion_constant(pss)
        errs.append(abs(info['c_from_growth'] / c - 1.0))
    assert errs[0] < 0.03, 'coarsest grid off by %.3f' % errs[0]
    for a, b in zip(errs, errs[1:]):
        assert b < 0.7 * a, \
            'the disagreement %.4f -> %.4f is not first-order; a ' \
            'difference that does not converge away is a defect, not ' \
            'discretisation' % (a, b)


def test_the_bordered_kronecker_solve_matches_its_closed_form():
    """`d = (vᵀK₁v)/(v·u)²` — the `O(n²)` contraction behind the `O(n⁴)` solve.

    Left-multiplying the bordered system by `(v⊗v)ᵀ` annihilates the
    singular block, so the growth rate never needed the Kronecker at all.
    They are the same quantity by construction, which makes a disagreement
    diagnostic rather than a precision question: it would mean the border
    pair is not the null pair.

    ⚠ AND THE DEFLATION IS WHAT MAKES EITHER COMPUTABLE. `I − M⊗M` has
    `σ_min = 2.3e-11` against a next singular value of 0.997 — a cleanly
    one-dimensional null space, which is why bordering with a single pair
    is the right repair. Bordered, `σ_min` comes back to 5.4e-02: nine
    orders recovered, the same shape as A7's deflated PAC solve.
    """
    _pss, _pac, _K, d, info = _osc_cov()
    assert info['d_residual'] < 1e-9, \
        'the solve and the closed form differ by %.3e; the border pair ' \
        'is not the null pair' % info['d_residual']
    assert info['sigma_min'] < 1e-8, \
        'sigma_min = %.3e; I - M kron M is supposed to be SINGULAR here ' \
        '-- if it is not, this circuit is not autonomous' % info['sigma_min']
    assert info['sigma_min_bordered'] > 1e-3, \
        'bordered sigma_min = %.3e; the deflation did not recover the '\
        'conditioning' % info['sigma_min_bordered']
    assert info['null_residual'] < 1e-8


def test_the_pair_inner_product_is_not_one_and_d_scales_with_the_tangent():
    """⚠ THE 2.31× THAT WAS CHASED AS A CODE DEFECT, PINNED AS A NUMBER.

    `ppv()` normalises on the FIRST BLOCK, `v[:m]·ẋ = 1` — right for a
    perturbation entering the first block, which is where an injected
    current lands and what every shipped path does. The FULL PAIR
    contraction is a different number, ≈0.663, so `(v·u)⁻² ≈ 2.27`.
    Assuming the pair product is 1 is exactly the error that produced a
    2.31× discrepancy between two Monte Carlo routes.

    ⚠ AND `d` ALONE IS NOT AN INVARIANT. Rescaling `u → s·u` sends
    `d → d/s²`; only the product `d·u uᵀ` is a property of the circuit.
    Asserted directly, because a future change to how the tangent is
    scaled would silently move `d` while leaving every structural check
    passing — `info['growth']` is what downstream code should read.
    """
    _pss, _pac, _K, d, info = _osc_cov()
    vu = info['pair_inner']
    assert 0.5 < abs(vu) < 0.9, \
        'v.u = %.4f; if this has become 1.0 the pair normalisation ' \
        'changed and every d is off by %.3f' % (vu, 1.0 / vu ** 2)
    u = info['tangent_pair']
    assert np.allclose(info['growth'], d * np.outer(u, u))
    ## the invariant, stated as a rescaling that must not move it
    s = 3.0
    assert abs((d / s ** 2) * np.outer(s * u, s * u)
               - d * np.outer(u, u)).max() < 1e-18


def test_the_oscillator_covariance_refuses_a_driven_circuit():
    """The mirror of `covariance`'s refusal, and it is not symmetry for its
    own sake. A driven circuit's `I − M⊗M` is nonsingular, so there is no
    null direction to border with and no walk to split off; bordering it
    anyway would return a `d` near zero and an arbitrary `K_orb`, which is
    a plausible wrong answer rather than an error.
    """
    import warnings
    cir = _rc_noisy()
    pss = PSS(cir, method='gear')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 400, refnode=gnd)
    with pytest.raises(ValueError, match='GROWS'):
        PAC(cir, toolkit=circuit.numeric).oscillator_covariance(pss)


def test_diffusion_constant_sees_noise_on_an_algebraic_row():
    """⚠ `diffusion_constant` used to return EXACTLY ZERO for an oscillator
    whose only noise is its series tank loss.

    The same physical loss, two equivalent representations, one answer --
    and cross-checked against `oscillator_covariance`, which reaches `CY`
    through the Lyapunov recursion instead of the PPV and was therefore
    right all along. The two routes disagreeing by the WHOLE quantity is
    how this was found.
    """
    cs, ps, pacs, _rs = _loss_osc('series', Q=8.0, npts=480)
    _cp, pp, pacp, _rp = _loss_osc('parallel', Q=8.0, npts=480)

    c_series = pacs.diffusion_constant(ps)
    c_parallel = pacp.diffusion_constant(pp)
    assert c_series > 0.0, \
        'noise on an algebraic row must not vanish; got %r' % c_series
    assert abs(c_series / c_parallel - 1.0) < 5e-3, \
        'series %.9e vs parallel %.9e' % (c_series, c_parallel)

    _K, _info = pacs.oscillator_covariance(ps)
    d = _info['d']
    dT = d / float(ps.period)
    assert abs(c_series / dT - 1.0) < 5e-3, \
        'the PPV route and the Lyapunov route must now agree: c %.9e, ' \
        'd/T %.9e' % (c_series, dT)


def _ghanta_tank(Q=16.0, cc=1.0, ll=1.0, psd=1e-6, npts=480):
    """An LC tank matching Ghanta, Li & Roychowdhury 2004 Lemma 5.2's premises.

    ODD-symmetric `i-v` (no even term) and a near-sinusoidal orbit, which the
    lemma requires.  `mu` is scaled with `C*w0` so `Q` means the same thing
    as `L` and `C` move.

    ⚠⚠ THE PREMISE GUARD WAS `rms/peak` AND IT MEASURED THE INTEGRATOR, NOT
    THE ORBIT (docs session, 2026-09-09, reproduced here).  0.7078470 at
    480 points is BDF-2 truncation: the true orbit's rms/peak is 0.7071072
    (DOP853 at 1e-13), a deviation of 2.8e-6 against the 7.4e-4 the guard
    read -- 265x the physical value, converging at exactly second order
    (0.70858 / 0.70785 / 0.70748 / 0.70729 at 240 / 480 / 960 / 1920).
    And the metric is SECOND order in the thing it guards: van der Pol's
    third harmonic is in quadrature with the fundamental, so it moves the
    peak only as `h3^2`, and an orbit with h3/h1 = 4.35 % (c 4.3 % off the
    lemma) passes a 2e-3 rms/peak gate that the assertion below holds to
    0.03 %.  The guard is now `h3/h1` from an rfft of the orbit's own
    samples: first order in the premise, grid-independent (1.243e-3 at
    480 and 960 points alike), no reference.  ⚠ The stored waveform has
    `steps + 1` samples (the endpoint repeats t = 0); the FFT takes the
    first `len(fp.steps)`.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    w0 = 1.0 / np.sqrt(ll * cc)
    mu = cc * w0 / (2 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=cc)
    cir['L'] = L('v', gnd, L=ll)
    cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    T = 2 * np.pi / w0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=250)
    assert pss.converged, 'C=%r L=%r did not converge' % (cc, ll)
    ## ⚠ waveform is (times, states): [0] is the TIME vector, [1] the states
    X = np.asarray(pss.waveform[1], dtype=float)
    row = X[0 if 0 < pss.irefnode else 1]
    amp = float(np.max(np.abs(row)))
    ns = len(pss.factored_period().steps)
    spec = np.abs(np.fft.rfft(row[:ns]))
    h3_h1 = float(spec[3] / spec[1])
    ## ⚠⚠ `psd / 2` TWICE, AND THE SECOND ONE IS THE POINT.  `noisePSD` and
    ## `_cy_reduced` are ONE-SIDED (a resistor's `4kT/R`, measured against
    ## the analytic value to every printed digit), while Ghanta's `N^2` is a
    ## TWO-SIDED density -- Winkler, Oberwolfach Report 18/2006 p.1160 gives
    ## Nyquist as `I_th = sqrt(2kT/R) xi(t)`, which is `2kT/R` two-sided and
    ## is exactly what `diffusion_constant` contracts (`cy/2`).  So the
    ## one-sided `psd` must be halved BEFORE entering the lemma.
    ##
    ## ⚠ Feeding the one-sided value gave a ratio that was CONSTANT at
    ## 0.49984 across 100x in `C` and 100x in `L` -- a perfect-looking gate
    ## that would have preserved the error forever while passing.  Halving
    ## it gives 0.99968 with no free parameter.  This campaign has already
    ## lost time to one factor of two that only `kT/C` could see, which is
    ## why the conversion is named here rather than absorbed.
    n_sq_two_sided = psd / 2.0
    lemma = (n_sq_two_sided / 2.0) * (ll / cc) / amp ** 2
    return cir, pss, PAC(cir, toolkit=circuit.numeric), amp, h3_h1, lemma


def test_the_lyapunov_route_matches_an_analytic_external_oracle():
    """⚠ THE FIRST FULLY EXTERNAL ORACLE FOR THE PHASE DIFFUSION CONSTANT.

    Ghanta, Li & Roychowdhury 2004 ASP-DAC Lemma 5.2, for an LC oscillator
    with an odd-symmetric `i-v` and a sinusoidal steady state:

        c = (N^2 / 2) (L / C) / A^2

    Every other check this codebase has on `c` is internal or shares the
    monodromy. This one shares nothing.

    ⚠ THE SWEEP IS THE TEST, NOT THE CONSTANT. A single operating point
    fixes only a scale factor, and a scale factor is exactly what a
    one-sided/two-sided PSD convention looks like. Sweeping `L` and `C`
    INDEPENDENTLY tests the functional form: `L/C` over four decades, at a
    constant ratio.

    ⚠ THE RATIO IS 1, NOT 1/2, ONCE THE PSD CONVENTION IS MATCHED. An
    earlier version of this gate fed the ONE-SIDED `noisePSD` into a lemma
    whose `N^2` is TWO-SIDED and asserted "constant at ~0.5" -- which passed,
    across four decades of `L/C`, while quietly carrying a factor of two.
    Winkler (Oberwolfach Report 18/2006 p.1160) gives Nyquist as
    `I_th = sqrt(2kT/R) xi(t)`, i.e. `2kT/R` TWO-SIDED, which is exactly the
    `cy/2` that `diffusion_constant` contracts. With the conversion made the
    ratio is 0.99968 with no free parameter, and the assertion is that it is
    ONE.
    """
    ratios = []
    for cc, ll in ((0.1, 1.0), (1.0, 1.0), (10.0, 1.0), (1.0, 0.1), (1.0, 10.0)):
        _cir, pss, pac, _amp, h3_h1, lemma = _ghanta_tank(cc=cc, ll=ll)
        ## ⚠ h3/h1, not rms/peak: see `_ghanta_tank`.  1.243e-3 here; 5e-3
        ## holds `c` to well under the 0.03 % the assertion below is at.
        assert h3_h1 < 5e-3, \
            'the lemma presumes a SINUSOIDAL orbit; h3/h1 came out %.3e' % h3_h1
        _K, _info = pac.oscillator_covariance(pss)
        d = _info['d']
        ratios.append((d / float(pss.period)) / lemma)
    lo, hi = min(ratios), max(ratios)
    assert abs(hi / lo - 1.0) < 1e-3, \
        'the ratio to the analytic oracle must be CONSTANT across L and C; ' \
        'it ranged %.6f to %.6f' % (lo, hi)
    assert abs(lo - 1.0) < 5e-3, \
        'and equal to ONE once the two-sided PSD convention is matched ' \
        '(Winkler, Oberwolfach 18/2006 p.1160); got %.6f' % lo


def test_diffusion_constant_should_not_depend_on_the_capacitance_scale():
    """⚠ `diffusion_constant` WAS WRONG BY EXACTLY `C^2`, and this is the sweep
    that found it -- against `oscillator_covariance`, which reaches `CY`
    through the Lyapunov recursion and never touches the PPV.

    The cause was contracting `CY`, an EQUATION-ROW covariance, against
    `samples` (`C^T v_1`) instead of `samples_eq` (`v_1`). Nothing caught it
    because every other fixture in this file uses `C = 1 F`, where the
    factor is exactly 1 -- section D shape 0i.
    """
    out = []
    for cc in (0.1, 1.0, 10.0):
        _cir, pss, pac, _amp, _rp, _lemma = _ghanta_tank(cc=cc, ll=1.0)
        _K, _info = pac.oscillator_covariance(pss)
        d = _info['d']
        out.append(pac.diffusion_constant(pss) / (d / float(pss.period)))
    lo, hi = min(out), max(out)
    assert abs(hi / lo - 1.0) < 1e-2, \
        'the PPV route and the Lyapunov route must agree at every C; the ' \
        'ratio ran %.6f to %.6f (a factor of %.1f, i.e. C^2)' \
        % (lo, hi, hi / lo)


def _tank_with_rc_probe(rpar, cpar, npts=480):
    """The lossy tank with a weakly-coupled noisy RC branch hung off node `x`.

    `R_par` is large against the tank's impedance so the branch does not load
    the oscillator, and `tau = R_par*C_par` is chosen rather than inherited --
    which is the whole point, because what follows is a statement about
    `tau/h`, not about `C`.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    Q = 8.0
    mu = 1.0 / (2 * np.pi * Q)
    rs = 0.2 * mu
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir.add_node('x')
    cir['L'] = L('v', 'x', L=1.0)
    cir['Rs'] = R('x', gnd, r=rs)
    cir.add_node('y')
    cir['Rpar'] = R('x', 'y', r=rpar)
    cir['Cpar'] = C('y', gnd, c=cpar)
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    T0 = 2 * np.pi
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / npts, x0=x0, maxiterations=250)
    assert pss.converged, 'rpar=%r cpar=%r' % (rpar, cpar)
    names = [str(nd) for nd in cir.nodes]
    iy = names.index('y')
    iy = iy if iy < pss.irefnode else iy - 1
    K, _i = PAC(cir, toolkit=circuit.numeric).oscillator_covariance(pss)
    return float(np.asarray(K, dtype=float)[iy, iy]), T0 / npts


def test_the_orbital_covariance_reaches_kTC_when_the_mode_is_RESOLVED():
    """⚠ `K_orb` DOES carry `kT/C` — what it cannot carry is an UNRESOLVED mode.

    This started as a suspected third defect: adding a parasitic capacitor at
    the algebraic node left `K_orb` flat over three decades of `C_par` and a
    factor ~1e6 BELOW `kT/C_par`. The explanation is not a missing term. That
    node's time constant was `rs*C_par ~ 4e-9 s` against a timestep of
    `0.013 s` -- **six orders faster than the grid**, and a mode the
    discretisation cannot represent cannot reach its equilibrium.

    ⚠⚠ THE DISCRIMINATOR IS THAT THE RATIO DEPENDS ON `tau/h` AND NOT ON `C`.
    At fixed `tau/h` it is identical across three decades of `C_par` -- so
    this is a resolution statement, not a scaling defect. Measured:

        tau/h = 152.79  ->  0.995110
        tau/h =  15.28  ->  0.953586
        tau/h =   1.53  ->  0.688971
        tau/h =   0.15  ->  0.214472
        tau/h =   0.02  ->  0.029182

    Monotone, and tending to `tau/h` itself once the mode is well below the
    grid. `kT/C` is an EXTERNAL anchor -- the same one that settled the
    `CY/2` convention -- so the top of that table is a real gate.
    """
    from pycircuit.circuit.constants import kboltzmann
    kT = kboltzmann * 300.0
    ## resolved: tau/h ~ 153
    kyy, h = _tank_with_rc_probe(rpar=1e5, cpar=2e-5)
    tau_over_h = (1e5 * 2e-5) / h
    assert tau_over_h > 100.0, 'this case must be well resolved; got %.1f' % tau_over_h
    assert abs(kyy / (kT / 2e-5) - 1.0) < 1e-2, \
        'a RESOLVED parasitic mode must reach kT/C: got %.6e against %.6e' \
        % (kyy, kT / 2e-5)

    ## the same tau/h at a different C -- the ratio must not move
    kyy2, _h = _tank_with_rc_probe(rpar=1e6, cpar=2e-6)
    r1 = kyy / (kT / 2e-5)
    r2 = kyy2 / (kT / 2e-6)
    assert abs(r2 / r1 - 1.0) < 1e-3, \
        'at fixed tau/h the ratio must be independent of C -- that is what ' \
        'makes this a RESOLUTION statement; got %.6f vs %.6f' % (r1, r2)

    ## unresolved: tau/h ~ 0.15, and the equilibrium must be largely absent
    kyy3, _h = _tank_with_rc_probe(rpar=1e2, cpar=2e-5)
    assert kyy3 / (kT / 2e-5) < 0.4, \
        'an UNRESOLVED mode must NOT reach kT/C, or this test shows nothing; ' \
        'got ratio %.6f' % (kyy3 / (kT / 2e-5))


def test_trap_oscillator_covariance_goes_through_the_twin_default_radau():
    """trap hands `oscillator_covariance` to the monodromy twin -- Radau by
    default since 2026-09-24 (TR-BDF2 before), Gear-2 selectable -- and the
    twin re-converges the SAME discrete orbit, so trap+twin matches a direct
    solve of the twin's method.

    ⚠ The NATIVE path refused until 2026-09-24 (`oscillator_covariance`
    bordered with `ppv()`'s width-m vectors, the trap pair map `2m x 2m`); it
    now runs on the pair's own null vectors -- first order, trap's own map's
    limit (`test_oscillator_covariance_runs_on_the_trapezoidal_pair_map`).
    The twin is what makes the default path accurate.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def build():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    T = 6.6634

    def solve(method, x0, mono=None):
        cir = build()
        pss = PSS(cir, method=method, reltol=1e-12)
        if mono is not None:
            pss.monodromy = mono
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 400, x0=np.asarray(x0),
                      maxiterations=100)
        assert pss.converged
        return pss, np.asarray(PAC(cir).oscillator_covariance(pss)[0],
                               dtype=float)

    ## default: trap hands off to the Radau twin.  Compared to a direct
    ## Radau solve SEEDED FROM TRAP'S OWN x0 -- the covariance at t=0 is a
    ## point on the orbit, so the two must be at the SAME PHASE to compare
    ## (seeding both elsewhere differs by O(h^2) of phase, ~1e-3 here, which
    ## is the orbit's covariance variation, not an error).
    tp, Kt = solve('trap', np.array([2.0, 0.0]))
    x0t = np.asarray(tp._period_state[1], dtype=float)
    _pd, Kd = solve('radau', x0t)
    assert np.linalg.norm(Kt - Kd) / np.linalg.norm(Kd) < 1e-6, \
        'the default trap twin and a direct Radau solve from the same ' \
        'seed differ by %.2e; the twin must re-converge the same orbit and ' \
        'injection' % (np.linalg.norm(Kt - Kd) / np.linalg.norm(Kd))
    ## the twin exists and is Radau by default
    assert tp.monodromy_twin() is not tp
    assert tp.monodromy_twin().par.method == 'radau'

    ## gear is selectable and matches a direct Gear-2 solve from the same seed
    _pg, Kg = solve('trap', np.array([2.0, 0.0]), mono='gear')
    x0g = np.asarray(_pg._period_state[1], dtype=float)
    _pdg, Kdg = solve('gear', x0g)
    assert np.linalg.norm(Kg - Kdg) / np.linalg.norm(Kdg) < 1e-6, \
        'monodromy=gear should match a direct Gear-2 solve from the same seed'

    ## native runs on the pair map's own null vectors (2m wide)
    tp.monodromy = 'native'
    _Kn, info_n = PAC(tp.cir).oscillator_covariance(tp)
    dn = info_n['d']
    assert np.shape(info_n['ppv_pair'])[0] == 2 * (tp.cir.n - 1)
    assert np.isfinite(dn) and dn > 0.0


def test_oscillator_covariance_takes_a_coloured_source_on_a_staged_oscillator():
    """Coloured `oscillator_covariance` on a STAGED oscillator (2026-09-26;
    Andreas: "Do 1").  The transverse responses are the total map's
    (`EventColumns.total_matrix`), with the source's own motion of the
    crossings in the right-hand side and the crossings' motion at FIXED
    time on the nodes -- `_forced_responses`' staged-oscillator path,
    replayed from the bounded part and projected.

    The gate: on the comparator relaxation oscillator the transverse path
    IS the projection of that verified full response, node by node -- 5e-13
    at 0.3 and 2.7 f0 (the pole's part it never forms projects to round-off
    there).  D2 exactly: a 1/f source added leaves `d` and `K_orb`
    bit-identical.

    ⚠ NOT gated by the white limit on this fixture, and measured why: a
    white source through the band route against the Lyapunov route reads
    -0.9 / -0.6 / -0.25 / -0.18 % at t = 0 (100 / 200 / 400 / 800 points)
    but 18 / 9.6 / 6.1 / 3.8 % at worst in the switch's ON phase, where the
    capacitor discharges through 10 ohm with `tau_on` about one step: the
    band a WHITE source has above the grid's Nyquist, (2/pi) f_on / f_N
    (7 / 3.5 % predicted at 400 / 800), is what the band route leaves out.
    The cycle mean converges more slowly (3.2 / 2.0 / 1.8 / 1.3 %); the
    Lyapunov route is itself first order on switched circuits (see
    `covariance`).  A coloured source has no such band.  (The Lyapunov
    route's own projection does not leak: `|Pi u_j|/|u_j|` ~ 1e-9.)"""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    res = {}
    for flick in (False, True):
        osc = _comparator_relaxation_oscillator()
        if flick:
            osc['fl'] = _Flicker('c', gnd, i=0.0, noisePSD=1e-22, fref=1.0)
        seed, To = _relaxation_oscillator_seed(osc)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss = PSS(osc, method='radau', reltol=1e-9)
            pss.solve(period=To, timestep=To / 200, x0=seed, maxiterations=100,
                      state_events=True)
        assert pss.converged and getattr(pss, '_event_columns', None) is not None
        pac = PAC(osc, toolkit=circuit.numeric)
        f0 = 1.0 / float(pss.period)
        if flick:
            with pytest.warns(RuntimeWarning, match="WHITE sources' alone"):
                res[flick] = pac.oscillator_covariance(pss, samples=True,
                                                       colour_fmin=1e-4 * f0)
        else:
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                res[flick] = pac.oscillator_covariance(pss, samples=True)
            ## the transverse path against the projected full response
            names = [str(x) for x in osc.nodes if str(x) != 'gnd!']
            ic = names.index('c')
            host = pss._lyapunov_host()
            fp = host._state_map()
            e = np.zeros(osc.n - 1, dtype=complex)
            e[ic] = 1.0
            Pi = pac._node_projectors(host)
            for nu in (0.3 * f0, 2.7 * f0):
                ys, _d = pac._forced_responses(host, fp, np.array([nu]), e)
                full = np.einsum('jab,jb->ja', Pi[:len(ys[0])], np.asarray(ys[0]))
                try:
                    tv, _n = pac._transverse_responses(host, fp, np.array([nu]), e)
                finally:
                    pac._transverse_cache = None
                n = min(len(full), len(tv[0]))
                err = np.max(np.abs(np.asarray(tv[0])[:n] - full[:n])) \
                    / np.max(np.abs(full[:n]))
                assert err < 1e-10, (nu / f0, err)
    (Kw, _iw), (Kc, ic_) = res[False], res[True]
    dw, dc = _iw['d'], ic_['d']
    assert dc == dw and np.array_equal(Kc, Kw)
    kc = np.asarray(ic_['coloured_samples'])
    assert np.all(np.diagonal(kc, axis1=1, axis2=2) >= -1e-30)
    assert np.max(kc) > 0.0


def test_the_transverse_band_integral_is_the_lyapunov_route_for_a_white_source():
    """The coloured oscillator path's WHITE LIMIT (2026-09-25): a white source
    pushed through the frequency-domain transverse route (as a component of
    exponent 0) against the Lyapunov route's ``Pi K_orb Pi^T``, node by
    node, on the same radau van der Pol.  The frequency route stops at the
    grid's Nyquist, and a WHITE source's transverse response falls as
    1/nu^2, so its missing tail is FIRST order in h: -5.2e-4 / -2.6e-4 /
    -1.3e-4 at 100 / 200 / 400 points, the Lyapunov route flat to seven
    digits.  Richardson in N removes it: +5e-6 (100/200), +2e-6 (200/400).
    (A 1/f source's tail is h^2, a Lorentzian's h^3 -- and white sources
    never take this route outside this check.)"""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    Q = 8.0
    mu = 1.0 / (2.0 * np.pi * Q)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    rel = {}
    for npts in (100, 200):
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0))
        c['n'] = IS('v', gnd, i=0.0, noisePSD=1e-8)
        pss = PSS(c, method='radau', reltol=1e-12)
        x0 = np.zeros(c.n - 1)
        x0[0] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=300)
            pac = PAC(c, toolkit=circuit.numeric)
            _K, info = pac.oscillator_covariance(pss, samples=True)
            host = pss._lyapunov_host()
            fp = host._state_map()
            counts, states = pac._injection_points(host, fp)
            f0 = 1.0 / float(host.period)
            W = psd_sqrt(pac._noise_components(host, states).cy_at_states(
                2 * np.pi * f0))
            col = {'fp': fp, 'counts': counts, 'comps': [(('n',), W, 0.0)],
                   'w1': 2 * np.pi * f0, 'perband': [], 'state0': states[0],
                   'fmin': 1e-7 * f0, 'fmax': 0.5 * len(fp.steps) * f0,
                   'ppd': 40}
            try:
                Kf, _ = pac._coloured_covariance(
                    host, col, c.n - 1, c.n - 1,
                    responses=pac._transverse_responses,
                    lines=pac._orbital_lines(host, col['fmin'], col['fmax']))
            finally:
                pac._transverse_cache = None
        tr = np.asarray(info['transverse_samples'])
        n = min(len(tr), len(Kf))
        rel[npts] = float(np.mean(Kf[:n, 0, 0]) / np.mean(tr[:n, 0, 0]) - 1.0)
    assert -1e-3 < rel[100] < -1e-4 and -5e-4 < rel[200] < -5e-5, rel
    assert abs(2.0 * rel[200] - rel[100]) < 2e-5, rel


def _a11_osc_chain(psd_tank=1e-6, psd_buf=1e-6, nstage=3, cb=0.5):
    """A11's AUTONOMOUS fixture: a van der Pol tank driving tanh buffers.

    `BSource` reads the tank voltage without drawing current, so the buffers
    do not load it and the orbit is van der Pol's own -- the autonomous
    analogue of `_a8_buffer_chain`.

    ⚠ `cb = 0.5` IS NOT A FREE CHOICE.  At `cb = 0.1` (T/tau = 66.6) the
    FROZEN free-period Jacobian goes singular: the fast buffer rows of
    `I - M` go near-degenerate and the frozen pin lands there.
    `phase_rule='reselect'` rescues it onto the SAME orbit (periods agree to
    4.44e-16), but `cb = 0.5` converges on the shipped defaults and needs no
    opt-in flag.  ⚠ The library's error blames "a seed BELOW the
    fundamental"; that is NOT this cause -- the seed is the measured period
    and the bare tank converges from it.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    if psd_tank > 0:
        c['n'] = IS('v', gnd, i=0.0, noisePSD=psd_tank)
    prev = 'v'
    for k in range(nstage):
        nd = 'o%d' % k
        c.add_node(nd)
        c['B%d' % k] = BSource(prev, gnd, nd, gnd,
                               i_func=lambda u: 1.0 * np.tanh(2.0 * u))
        c['R%d' % k] = R(nd, gnd, r=1.0, noisy=False)
        c['C%d' % k] = C(nd, gnd, c=cb)
        if psd_buf > 0:
            c['nb%d' % k] = IS(nd, gnd, i=0.0, noisePSD=psd_buf)
        prev = nd
    return c


def _a11_solved(psd_tank=1e-6, psd_buf=1e-6, npts=240, method='gear'):
    """`(cir, pss, reduced index of the last buffer, crossing time)`."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _a11_osc_chain(psd_tank, psd_buf)
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.664052486, timestep=6.664052486 / npts, x0=x0,
                  maxiterations=80)
    assert pss.converged, 'the A11 autonomous fixture did not converge'
    ## ⚠ this circuit has an INDUCTOR, so the reduced index is NOT the
    ## position among non-gnd nodes; it is the inverse of the lifting rule
    full = cir.get_node_index('o2')
    red = full if full < cir.get_node_index(gnd) else full - 1
    Xw = np.asarray(pss.waveform[1], dtype=float)
    grid = np.asarray(pss.factored_period().times, dtype=float)[:Xw.shape[1]]
    v = Xw[full]
    mid = 0.5 * (v.max() + v.min())
    j = [k for k in range(2, len(v) - 2)
         if (v[k - 1] - mid) < 0 <= (v[k] - mid)][0]
    tc = grid[j - 1] + (mid - v[j - 1]) / (v[j] - v[j - 1]) * (grid[j] - grid[j - 1])
    return cir, pss, red, tc, v


def test_the_additive_edge_jitter_of_an_oscillator_is_the_projected_bounded_covariance():
    """A11's second half: the NON-ACCUMULATING edge jitter of a free-running
    oscillator -- the number `c` does not contain.

    `sampled_noise` refuses an autonomous PSS (a diffusing phase has no
    sampling instant), so `jitter_metrics` cannot be used here.  The split
    `oscillator_covariance` returns supplies it instead: `n d u_j u_jᵀ` is the
    walk, `P(t_j)` is bounded, and the additive jitter is the bounded part
    OBLIQUELY PROJECTED, `Π P Πᵀ` with `Π = I − u_j v_jᵀ/(v_jᵀu_j)`.

    ⚠⚠ THE SCOPE THIS CLOSES WAS MIS-STATED BY ME.  A11 said the remaining
    half was "bias-dependent, non-stationary delay modulation, outside the
    modulated-stationary support".  Andreas concurred it was wrong
    (2026-09-17): A8's own entry already had the LDO injected as a
    SUPPLY-NODE NOISE SOURCE in "one ordinary autonomous PSS ... no multirate
    anything".  There is no exotic support problem; the gap was only that no
    route existed for the additive part.

    MEASURED 2026-09-17, over EIGHTY Monte Carlo seed-runs -- noisy transients
    carrying no PSS, no adjoint and no Lyapunov solve:

        MC / analysis = 1.0066 ± 0.0102   (0.64σ from 1.000)

    over 124 seed-runs on three grids, once both sides use the same orbit.
    ⚠ SUPERSEDED 2026-09-29: that analysis was the TWO-sided `A`, whose
    k-lag intercept drops the transverse-phase cross term; a COMMITTED Monte
    Carlo (`benchmarks/oscillator_edge_jitter_probe.py`) excludes it at 8.9σ
    and agrees with the exact law now returned -- see
    `test_the_oscillator_edge_jitter_is_the_exact_law_the_monte_carlo_measures`.
    The 2026-09-17 harness was never committed, so its statistic cannot be
    reproduced.

    ⚠⚠ AN EARLIER VERSION OF THIS DOCSTRING CLAIMED A 4σ DEFECT AND IT WAS
    WRONG. An 80-seed campaign gave 1.0571 ± 0.0141 and I recorded it as "the
    route under-predicts by 5.7 %, grid-independent, therefore the analysis's
    problem". The fault was the COMPARISON: `σ_t = √var/slew`, so it goes as
    `1/slew²`, and the Monte Carlo crosses on an EULER orbit while the method
    divides by the PSS's GEAR slope. Those orbits differ at finite h — Euler's
    slope at the crossing is low by 2.717 / 1.361 / 0.686 % at npts
    240 / 480 / 960, halving as first order requires — so squared they predict
    1.0566 / 1.0278 / 1.0139. Dividing each grid by its own measured
    correction collapses the raw 1.0493 / 1.0513 / 1.0165 onto 0.9931 / 1.0229
    / 1.0025, pooling to 1.0066 ± 0.0102.

    ⚠ THE LESSON OUTLIVES THE NUMBER: "grid-independent" was asserted from two
    points 0.99σ apart. Two noisy points cannot distinguish FLAT from HALVING,
    and halving is what it was doing — an absence of evidence read as evidence
    of a property, with a structural conclusion built on top.

    WHAT THIS TEST GATES is the STRUCTURE and the invariants below, plus the
    exclusion of a FACTOR (a factor of 2 in the variance is >15σ away). It
    deliberately does not assert the ratio, which needs ~100 transient runs.
    """
    cir, pss, red, tc, _v = _a11_solved()
    pac = PAC(cir, toolkit=circuit.numeric)
    r = pac.oscillator_edge_jitter(pss, red, tc)

    ## 1. the shipped two-anchor loop must still close on this 5-state
    ##    fixture -- a forward Lyapunov recursion against an adjoint-replayed
    ##    PPV, sharing only the CY/2 convention
    c_ref = float(pac.diffusion_constant(pss))
    assert abs(r['c'] / c_ref - 1.0) < 5e-3, (r['c'], c_ref)

    ## 2. the projection is not decoration, and its size is a property of the
    ##    SOURCE MIX: the phase direction is what tank noise drives, so a
    ##    quiet tank makes it vanish.  Measured 0.3728 here against 0.00034
    ##    when the tank is 1000x quieter (the ONE-sided projection since
    ##    2026-09-29; the two-sided one read 0.1611 / 0.0002).
    assert 0.30 < r['projection_share'] < 0.45, r['projection_share']
    _c2, p2, red2, tc2, _v2 = _a11_solved(psd_tank=1e-9)
    r2 = PAC(_c2, toolkit=circuit.numeric).oscillator_edge_jitter(p2, red2, tc2)
    assert r2['projection_share'] < 0.01, \
        'a quiet tank should leave almost nothing for the projection to ' \
        'remove, got %.4f' % r2['projection_share']

    ## 3. ⚠ GRID INDEPENDENCE, and it is a real property only because the
    ##    slope is taken at the REQUESTED INSTANT.  Differentiating at the
    ##    nearest grid point instead reads 1.485472 on this grid against a
    ##    converged 1.528 -- 2.8 % low, 5.6 % in the variance -- and that
    ##    error masqueraded as the covariance converging.
    _c3, p3, red3, tc3, _v3 = _a11_solved(npts=480)
    r3 = PAC(_c3, toolkit=circuit.numeric).oscillator_edge_jitter(p3, red3, tc3)
    assert abs(r3['A'] / r['A'] - 1.0) < 1e-2, (r['A'], r3['A'])
    assert abs(r3['slew'] / r['slew'] - 1.0) < 5e-3, (r['slew'], r3['slew'])

    ## 4. the k-lag law: Var -> c k T + 2A, so k_cycle rises from just above
    ##    sqrt(2A) and the walk eventually dominates
    assert r['k_cycle'][0] >= np.sqrt(2.0 * r['A'])
    assert np.all(np.diff(r['k_cycle']) > 0)

    ## 5. refusal on a DRIVEN circuit.
    ## ⚠ THE FIXTURE IS BUILT OUTSIDE THE `raises` BLOCK, DELIBERATELY.  With
    ## the construction inside it, a ValueError from the FIXTURE would satisfy
    ## the context and the method under test would never be called at all --
    ## a vacuous pass, and a green dot cannot tell the two apart.  Only the
    ## call being tested belongs inside.
    cirq, pssq, _io, _pacq, _T = _sampler_fixture(
        lambda c: c.__setitem__('S0', _sw()))
    assert not getattr(pssq, 'autonomous', False), \
        'the driven-refusal control is not driven, so it proves nothing'
    with pytest.raises(ValueError, match='covariance that GROWS'):
        PAC(cirq, toolkit=circuit.numeric).oscillator_edge_jitter(pssq, 0, 0.0)


def test_the_oscillator_consumers_read_the_total_map_on_a_staged_solve():
    """Item 1 of the 2026-09-22 list: `floquet_modes` diagonalised the
    FIXED-GRID map of a staged solve and sampled its modes along the orbit
    without the event nodes' costate injections; `oscillator_covariance`
    bordered `I - M kron M` with that map while its `u`, `v` were the total
    map's.  Now the total map (`M + P_theta dtheta/dx_0`) and the
    injections (`_event_costate_injection`, factored out of `ppv`).
    Measured on the staged comparator oscillator (compact switch, 200
    points): the multipliers move 2e-6 onto the total's (1.998273e-2; the
    exact 1.998264e-2); oscillator_covariance's bordered residual 1.06e-7
    -> 5.1e-15; the second mode's `C^T q` along the orbit matches the exact
    saltation propagation's direction to 6 digits -- with or without the
    injections, since with the transition inside the landed window the
    fixed-grid map carries the switching and the injections are the
    per-mille sliver (as everywhere today).  Pinned: multipliers equal to
    eig(M_tot) to 1e-9, the bordered residual below 1e-11, the mode's
    direction above 0.99999 at nodes before and after the crossings."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    names = [str(n_) for n_ in cir.nodes]
    seed, Tl = _relaxation_oscillator_seed(cir)
    q = PSS(cir, method='radau', reltol=1e-9)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100, state_events=True)
    assert q.converged and q._event_columns is not None
    fp = q.factored_period()
    n = fp.width
    Md = np.column_stack([np.asarray(fp.matvec(e), dtype=float) for e in np.eye(n)])
    Mt = Md + np.asarray(q._event_columns['P_end'], dtype=float) @ np.asarray(q._event_sensitivity, dtype=float)
    lam_t = np.sort(np.abs(np.linalg.eigvals(Mt)))[::-1]
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        fm = q.floquet_modes(nmodes=2)
    lam_fm = np.sort(np.abs([mm['lam'] for mm in fm]))[::-1]
    assert np.max(np.abs(lam_fm - lam_t[:2]) / lam_t[:2]) < 1e-9, (lam_fm, lam_t[:2])
    pac = PAC(cir, toolkit=circuit.numeric)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _K, info = pac.oscillator_covariance(q)
    assert info['d_residual'] < 1e-11, info['d_residual']
    ## the exact second left mode along the orbit (saltation propagation)
    mdl = _exact_relaxation_oscillator_model()
    lam_ex, V_ex = np.linalg.eig(mdl.M.T)
    o = np.argsort(-np.abs(lam_ex))
    v2 = np.real(V_ex[:, o[1]])
    assert abs(lam_t[1] / abs(lam_ex[o[1]]) - 1.0) < 1e-4, (lam_t[1], lam_ex[o[1]])
    left2 = mdl.left_at(v2)
    tg = np.linspace(0.0, mdl.T, 20001)[:-1]
    orb = np.array([left2(x)[0] for x in tg])
    red = [nm for i, nm in enumerate(names) if i != q.irefnode]
    idx = [red.index(nm) for nm in ('c', 'fb0', 'fb1')]
    Xw = np.asarray(q.waveform[1], dtype=float)
    sel = [names.index(nm) for nm in ('c', 'fb0', 'fb1')]
    qs = np.asarray(fm[1]['q'])
    for j in (10, 150):
        xs = Xw[sel, j]
        k = int(np.argmin(np.linalg.norm(orb - xs, axis=1)))
        _x, ve = left2(float(tg[k]))
        Cj = np.asarray(q._C_at(np.delete(Xw[:, j], q.irefnode)), dtype=float)
        w = (Cj.T @ np.real(qs[:, j]))[idx]
        cos = abs(w @ ve) / (np.linalg.norm(w) * np.linalg.norm(ve))
        assert cos > 0.99999, (j, cos)


def test_oscillator_covariance_runs_on_the_trapezoidal_pair_map():
    """`oscillator_covariance` on trap's OWN map (`monodromy='native'`) was
    refused: the plain trapezoidal state is the pair `(x, iq)`, its map
    `2m x 2m`, and the bordered solve needs THAT map's null vectors where
    `ppv()` gives width `m`.  They follow from the state map's: the map
    re-seeds `iq` at every period start, so its last `m` columns are zero,
    and the null vectors are `[v; 0]` and `M[:, :m] u`.  Measured on van der
    Pol against radau at 800 points: d -2.5 % at 200 points (-1.1 % at 800)
    -- first order, trap's own map's limit, as euler's is -- where the
    trbdf2 twin reads -7.6e-5.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    d_ref = 2.5756770857716342e-08          # radau, 800 points
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    out = {}
    for mono in ('native', 'trbdf2'):
        cir = _vdp_asym()
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        p = PSS(cir, method='trap', reltol=1e-10)
        p.monodromy = mono
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]),
                    maxiterations=60)
            _K, info = PAC(cir, toolkit=circuit.numeric).oscillator_covariance(p)
            d = info['d']
        out[mono] = (d, info)
    d, info = out['native']
    assert np.shape(info['ppv_pair'])[0] == 4 and info['d_residual'] < 1e-4, info['d_residual']
    assert abs(d / d_ref - 1.0) < 0.05, d / d_ref - 1.0
    assert abs(out['trbdf2'][0] / d_ref - 1.0) < 1e-3, out['trbdf2'][0] / d_ref - 1.0


def test_the_oscillator_covariance_and_mode_weights_return_one_shape():
    """Item #15 of `doc/pac_noise_conventions.md` (2026-09-29):
    `oscillator_covariance` returns `(K_orb, info)` -- it returned
    `(K_orb, d, info)` -- with the growth rate `info['d']` and the bounded
    part per node `info['samples']` (was `'orbital_samples'`);
    `orbital_mode_weights` returns `(cw, info)` with `info['modes']` and
    `info['K']` (was `(cw, modes, K)`)."""
    _cir, pss, pac = _solve_vdp_noise(npts=240)
    res = pac.oscillator_covariance(pss, samples=True)
    assert len(res) == 2
    _K_orb, info = res
    assert 'orbital_samples' not in info
    assert {'d', 'samples', 'growth_samples', 'times',
            'transverse_samples'} <= set(info)
    assert info['d'] > 0.0 and len(info['samples']) == len(info['growth_samples'])
    w = pac.orbital_mode_weights(pss)
    assert len(w) == 2 and set(w[1]) == {'modes', 'K'}


def test_the_oscillator_edge_jitter_is_the_exact_law_the_monte_carlo_measures():
    """#17 B4 of `doc/pac_noise_conventions.md` (Andreas: "Exact k_cycle +
    fix A", 2026-09-29).  The exact linear k-lag law of the edge time,

        Var_k = e^T (2 P_j + k G_j - M_j^k P_j - P_j M_j^k^T) e / s^2,

    is `k_cycle`, and its large-k intercept `2 A` with the ONE-sided
    projection `A = e^T Pi P_j e / s^2`.  ⚠ Until 2026-09-29 the method
    returned `k_cycle_bound = sqrt(c k T + 2 A)` with the TWO-sided
    `A = e^T Pi P_j Pi^T e / s^2`, which drops the transverse-phase cross term
    (X/A = -0.150 here).  A committed Monte Carlo
    (`benchmarks/oscillator_edge_jitter_probe.py`: noisy radau transients at
    the PSS's own step, 8 seeds x 720 periods, 5760 crossings) measured
    Var_k = 1.5001e-6 / 2.0227e-6 / 3.058e-6 / 5.05e-6 (+- 2.2e-8 / 3.0e-8 /
    5.5e-8 / 2.0e-7) at k = 1 / 2 / 4 / 8: the exact law at -1.2 to -1.7 sigma,
    the old k = 1 value (1.6969e-6) at -8.9 sigma.  Radau, 240 points: the
    law's own numbers are pinned (the probe's), and each sits within 3 sigma
    of the Monte Carlo."""
    cir, pss, red, tc, _v = _a11_solved(method='radau')
    r = PAC(cir, toolkit=circuit.numeric).oscillator_edge_jitter(pss, red, tc)
    assert 'k_cycle_bound' not in r
    kc2 = r['k_cycle'] ** 2
    assert abs(kc2[0] / 1.52740e-6 - 1.0) < 1e-4, kc2[0]
    assert abs(2.0 * r['A'] / 9.86462e-7 - 1.0) < 1e-4, r['A']
    for k, mc, se in ((1, 1.5001e-6, 2.22e-8), (2, 2.0227e-6, 2.95e-8),
                      (4, 3.0584e-6, 5.53e-8), (8, 5.0497e-6, 1.99e-7)):
        assert abs(kc2[k - 1] - mc) < 3.0 * se, (k, kc2[k - 1], mc)
    ## the law tends to the walk plus the one-sided intercept
    T = float(pss.period)
    assert abs((kc2[0] - r['c'] * T) / (2.0 * r['A']) - 1.0) < 1e-2


def _vdp_colour_pair(kind, npts=200):
    """van der Pol (Q = 8, radau) with ONE physical noise two ways: a
    Lorentzian `IS(noiseTau = 0.3 T)` on the tank ('coloured'), or white
    noise through an RC into a linear `BSource` on it ('filtered') -- the
    exact realisation (one added state).  `(cir, pss, reduced index of the
    tank, its rising mid-level crossing)`."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    Rf, Cf, g, Pw = 1.0, 0.3 * T, 1e-2, 1e-4
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    if kind == 'coloured':
        c['n'] = IS('v', gnd, i=0.0, noisePSD=g * g * Pw * Rf * Rf,
                    noiseTau=Rf * Cf)
    else:
        c.add_node('f')
        c['nw'] = IS('f', gnd, i=0.0, noisePSD=Pw)
        c['rf'] = R('f', gnd, r=Rf)
        c['cf'] = C('f', gnd, c=Cf)
        c['gm'] = BSource('f', gnd, gnd, 'v', i_func=lambda u, _g=g: _g * u)
    pss = PSS(c, method='radau', reltol=1e-12)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=300)
    assert pss.converged
    full = c.get_node_index('v')
    red = full if full < c.get_node_index(gnd) else full - 1
    Xw = np.asarray(pss.waveform[1], float)
    grid = np.asarray(pss.factored_period().times, float)[:Xw.shape[1]]
    v = Xw[full]
    mid = 0.5 * (v.max() + v.min())
    j = next(k for k in range(2, len(v) - 2)
             if (v[k - 1] - mid) < 0 <= (v[k] - mid))
    tc = grid[j - 1] + (mid - v[j - 1]) / (v[j] - v[j - 1]) * (grid[j] - grid[j - 1])
    return c, pss, red, float(tc)


def test_the_oscillator_edge_jitter_takes_a_coloured_source_as_its_realisation_does():
    """#17 B2 (2026-09-29): a COLOURED source in `oscillator_edge_jitter`.
    Its transverse part joins `A`; its PHASE is the colour fold's
    INSTANT-SPECIFIC increment (`_lineshape.edge_increment`), not the
    stationary structure function -- a coloured source has memory, and at
    this edge the stationary form reads 1.003e-9 against the increment's
    7.53e-10 at k = 1 (400 points).
    Gated on the EDGE TIME, which is physical in any coordinates: the element
    against its exact white realisation (white noise through an RC into the
    tank) read by the exact law (B4, Monte-Carlo validated).  Predicted
    <= 0.5 %; measured +1.4e-3 / +2.4e-4 / ... / -4.3e-4 at k = 1..8 on 400
    points, within 3e-3 here (200 points, the band capped at 20 f0).
    ⚠ The realisation's own intercept is NEGATIVE (the RC state carries the
    noise the phase later takes up: anti-correlated), so it has no additive
    variance -- warned, `sigma_t` nan, `k_cycle` exact."""
    cw, pw, rw, tw = _vdp_colour_pair('filtered')
    with pytest.warns(RuntimeWarning, match='NEGATIVE'):
        w = PAC(cw, toolkit=circuit.numeric).oscillator_edge_jitter(pw, rw, tw)
    assert w['A'] < 0.0 and np.isnan(w['sigma_t']) and w['band'] is None
    ce, pe, re_, te = _vdp_colour_pair('coloured')
    f0 = 1.0 / float(pe.period)
    pac = PAC(ce, toolkit=circuit.numeric)
    with pytest.raises(TypeError, match='COLOURED'):
        pac.oscillator_edge_jitter(pe, re_, te)
    e = pac.oscillator_edge_jitter(pe, re_, te, colour_fmin=1e-6 * f0,
                                   colour_fmax=20.0 * f0, points_per_decade=20)
    assert e['A'] > 0.0 and e['projection_share'] is None
    assert np.all(e['coloured_phase_variance'] > 0.0)
    r = e['k_cycle'] ** 2 / w['k_cycle'] ** 2 - 1.0
    assert np.max(np.abs(r)) < 5e-3, r


def test_the_oscillator_edge_jitter_takes_a_modulated_flicker_as_its_realisation_does():
    """#17 B2, the ORBIT-MODULATED 1/f path (2026-09-29): a flicker current
    whose level follows the tank voltage (`k V_v flicker_noise(1)`, signed
    amplitudes: a POWER-LAW group of columns that move along the orbit)
    against its realisation -- a stationary 1/f source on an algebraic node
    times the tank voltage through a multiplier, no added state
    (`_orbit_modulated_vdp`, asymmetric at a = 0.3 so the 1/f up-converts).
    Two circuits, two routes into the same physics: the edge jitter's
    coloured transverse band integral and the fold's phase increment agree
    to rounding -- measured 3.7e-10 on `k_cycle^2` at every k, 3.7e-10 on
    `A`, 4.0e-10 on the increment (predicted <= 1e-6 / 1e-5).
    ⚠ The fixture's unit noise level suits the SPECTRAL ratios it was built
    for; as a timing it is sigma_t = 19 s, 2.8 periods, which the method's
    first-order guard refuses -- correctly.  Hence `level = 1e-8` here
    (sigma_t 1.9e-3 s)."""
    out = {}
    for kind in ('flicker', 'flicker_ref'):
        c, pss, pac, _ov = _orbit_modulated_vdp(kind, method='radau', a=0.3,
                                                level=1e-8, npts=200)
        f0 = 1.0 / float(pss.period)
        Xw = np.asarray(pss.waveform[1], float)
        grid = np.asarray(pss.factored_period().times, float)[:Xw.shape[1]]
        full = c.get_node_index('v')
        red = full if full < c.get_node_index(gnd) else full - 1
        v = Xw[full]
        mid = 0.5 * (v.max() + v.min())
        j = next(k for k in range(2, len(v) - 2)
                 if (v[k - 1] - mid) < 0 <= (v[k] - mid))
        tc = grid[j - 1] + (mid - v[j - 1]) / (v[j] - v[j - 1]) \
            * (grid[j] - grid[j - 1])
        out[kind] = pac.oscillator_edge_jitter(
            pss, red, tc, colour_fmin=1e-4 * f0, colour_fmax=20.0 * f0,
            points_per_decade=20)
    e, r = out['flicker'], out['flicker_ref']
    T = 1.0 / f0
    assert 0.0 < e['sigma_t'] < 1e-3 * T
    assert np.all(e['coloured_phase_variance'] > 0.0)
    assert np.max(np.abs(e['k_cycle'] ** 2 / r['k_cycle'] ** 2 - 1.0)) < 1e-8
    assert abs(e['A'] / r['A'] - 1.0) < 1e-8
    assert np.max(np.abs(e['coloured_phase_variance']
                         / r['coloured_phase_variance'] - 1.0)) < 1e-8
