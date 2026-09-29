"""Shooting tests: shooting covariance.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_elements import (_NuMult,
    _SwitchHdl)
from pycircuit.circuit.shooting._noise_components import separable
from pycircuit.circuit.tests._shooting_fixtures import _noise_seam
from pycircuit.circuit.tests._shooting_fixtures import (_Flicker,
    _KB,
    _ModLorentzCtl,
    _TEMP,
    _comparator_relaxation_oscillator,
    _jitter_sampler,
    _mixed_exponent_rc,
    _pow_flicker,
    _pwm_loop,
    _q20_rlc,
    _rc_noisy,
    _relaxation_oscillator_seed,
    _sampler_fixture_method,
    _solve_slow,
    _sw)


def _cov_ratio(npts, Cval=1e-7, per=1e-3):
    import warnings
    from pycircuit.circuit.constants import kboltzmann
    circuit.default_toolkit = circuit.numeric
    cir = _rc_noisy(Cval=Cval, per=per)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    K0 = pac.covariance(pss)
    irn = pss.irefnode
    k = cir.get_node_index(cir.get_node('b'))
    k = k - 1 if k > irn else k
    T = float(circuit.defaultepar.T)
    return K0[k, k] / (kboltzmann * T / Cval), K0


def test_the_periodic_covariance_converges_to_kTC():
    """⚠ AN EXACT, FAMOUS, INDEPENDENT ANSWER: `Var(v_C) = kT/C`.

    Independent of `R`, which makes it a strong check — nothing about the
    resistor, the drive or the grid should appear in it. The covariance
    here is built from per-step Lyapunov accumulation and one Kronecker
    solve, and shares no machinery with the closed form it is checked
    against.

    ⚠ THE ASSERTION IS THE RATE. Measured 0.931 / 0.964 / 0.982 / 0.991 at
    100 / 200 / 400 / 800 points — the error halving each doubling, which
    is the O(h) a piecewise-constant approximation to white noise gives.
    A wrong constant would sit at a fixed offset and never move; this is
    what distinguishes the two.

    ⚠ AND IT SETTLED A FACTOR OF TWO BY MEASUREMENT. With the full `CY` the
    same sequence converges to 2.0, not 1.0 — so `CY` is a ONE-SIDED
    density and the per-step injection carries `CY/2`. Argument alone
    would not have decided it; the two candidate conventions differ by
    exactly the factor the test resolves.
    """
    rs = [_cov_ratio(n)[0] for n in (100, 200, 400)]
    assert rs[-1] > 0.95, \
        'the covariance is %.4f of kT/C at 400 points; it should be ' \
        'approaching 1' % rs[-1]
    errs = [abs(1.0 - r) for r in rs]
    for a, b in zip(errs, errs[1:]):
        assert 1.6 < a / b < 2.6, \
            'the error falls %.2fx per doubling (%s), not the ~2x of O(h). ' \
            'A wrong CONSTANT would not move at all' % (a / b, errs)


def test_the_covariance_grid_must_resolve_the_noise_bandwidth():
    """⚠ A PRECONDITION, NOT AN ACCURACY NOTE — and it cost a wrong reading.

    The first attempt at the gate above came back at 0.517 and looked like
    a factor-of-two bug. It was not: the RC pole sat at 159 kHz while the
    grid's Nyquist was 100 kHz, so the DISCRETE system genuinely does not
    carry the noise the continuous one does. `kT/C` requires integrating
    past the pole.

    Reproduced here deliberately: the same circuit with the pole above
    Nyquist reads far low, and moving the pole below it recovers the
    answer. A `kT/C` that comes back low is the grid, not the code.
    """
    ## pole at 159 kHz, Nyquist at 100 kHz — under-resolved
    under, _K = _cov_ratio(200, Cval=1e-9)
    ## pole at 1.6 kHz, Nyquist at 100 kHz — resolved
    ok, _K2 = _cov_ratio(200, Cval=1e-7)
    assert under < 0.7, \
        'the under-resolved case reads %.4f; it is supposed to be visibly ' \
        'low, or this test is not demonstrating the precondition' % under
    assert ok > 0.95, 'the resolved case reads %.4f' % ok


def test_the_covariance_refuses_an_oscillator():
    """`I − M⊗M` is singular there, and that is the physics.

    The unit multiplier squares to one, the covariance grows without bound
    rather than settling, and that growth IS the phase diffusion. An
    oscillator's output noise is stationary, not cyclostationary — there is
    no object here to compute, so it refuses and names the right route.
    """
    _cir, pss = _solve_slow(None)
    with pytest.raises(ValueError, match='no periodic covariance'):
        PAC(pss.cir, toolkit=circuit.numeric).covariance(pss)


def test_the_covariance_samples_are_periodic_and_positive():
    """The time-varying statistic itself: `K(t)` over one period.

    Two structural properties that a wrong recursion breaks differently —
    it must return to itself after a period (that is what the Kronecker
    solve enforces), and every `K_j` must be a positive-semidefinite
    covariance, which nothing in the solve guarantees a priori.
    """
    _r, _K = _cov_ratio(100)
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _rc_noisy()
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
    K0, seq = PAC(cir, toolkit=circuit.numeric).covariance(pss, samples=True)

    rel = np.linalg.norm(seq[-1] - K0) / max(np.linalg.norm(K0), 1e-300)
    assert rel < 1e-8, \
        'K(T) differs from K(0) by %.3e -- the solve is supposed to ' \
        'enforce exactly that periodicity' % rel
    for j, K in enumerate(seq):
        w = np.linalg.eigvalsh(K)
        assert w.min() > -1e-12 * max(abs(w).max(), 1e-300), \
            'K at step %d has eigenvalue %.3e; a covariance cannot be ' \
            'negative definite' % (j, w.min())


def _driven_rlc_for_lyapunov(method, npts=200):
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    Lv, Cv, Rs = 1e-3, 1e-9, 10.0
    per = 2.0 * np.pi * np.sqrt(Lv * Cv)
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
    c['r'] = R('a', 'b', r=Rs)
    c['l'] = L('b', gnd, L=Lv)
    c['c1'] = C('b', gnd, c=Cv)
    c['n'] = IS('b', gnd, i=0.0, noisePSD=1e-18)
    import warnings
    pss = PSS(c, method=method, reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, x0=np.zeros(c.n - 1),
                  maxiterations=100, x0_unknown=False)
    assert pss.converged
    return c, pss


def test_plain_lyapunov_pieces_are_tied_to_the_shipped_monodromy():
    """The per-step maps of the PLAIN path, and the one gate that ties them
    to something already trusted: their product must BE `fp.matvec`.

    Euler's step state is `x` alone (`b = 0`), so the product is the
    monodromy directly. Trapezoidal's is the pair `(x, iq)`; the product
    applied to `(x, 0)` — the companion re-seeded, as the solve does —
    and read out on `x` must equal `fp.matvec` too. Gear is the control.

    ⚠ THE UN-RESET TRAP PAIR IS SINGULAR, measured: carrying `iq` across
    the period boundary puts a marginal `(−1)ⁿ` mode in `I − M⊗M`
    (`LinAlgError`), the obstruction this file records for every
    formulation that keeps the companion across a period. The period map
    re-seeds it, and that is what the tie below checks.
    """
    for method in ('euler', 'trap', 'gear'):
        cir, pss = _driven_rlc_for_lyapunov(method)
        pac = PAC(cir)
        fp = pss.factored_period()
        m = cir.n - 1
        As, Qs, K1, M, mm, n = pac._lyapunov_pieces(pss, 'covariance')
        Mfp = np.column_stack([np.asarray(fp.matvec(e), float)
                               for e in np.eye(fp.width)])
        P = np.eye(n)
        for A in As:
            P = A @ P
        got = P if n == fp.width else P[:m, :m]
        tie = float(np.linalg.norm(got - Mfp)) / float(np.linalg.norm(Mfp))
        assert tie < 1e-10, \
            '%s: the product of the per-step maps differs from fp.matvec ' \
            'by %.3e -- the Lyapunov pieces describe a different map from ' \
            'the one the solve converged on' % (method, tie)
        assert np.all(np.isfinite(K1)) and np.trace(K1) > 0.0, \
            '%s: the one-period noise accumulation is not positive' % method
        ## and the driven Kronecker solve must be NON-singular
        K0 = pac.covariance(pss)
        assert np.all(np.isfinite(K0))


def test_plain_covariance_reaches_kTC_like_gear_does():
    """`covariance` under euler and trap, against `kT/C` — a reference no
    integrator can influence.

    ⚠ THE FIXTURE MUST NOT DOUBLE-COUNT, and the first version did:
    adding an explicit `IS(noisePSD=4kT/R)` beside a resistor that already
    carries thermal noise made EVERY method read 2.0 kT/C (euler 1.939,
    gear 1.911 at 1600 points) — predicted before it was read, as the
    signature of double-counting rather than of any method, and it was.
    With the resistor's own noise only, every method converges on 1.0.
    """
    import warnings
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    kT = 1.380649e-23 * 300.0
    Rv, Cc = 1e3, 1e-9
    per = 100.0 * Rv * Cc

    def rc():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1e-3, freq=1.0 / per)
        c['r'] = R('a', 'b', r=Rv)
        c['c1'] = C('b', gnd, c=Cc)
        return c

    for method in ('euler', 'trap'):
        ratios = []
        for npts in (400, 1600):
            cir = rc()
            pss = PSS(cir, method=method, reltol=1e-11)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=per, timestep=per / npts,
                          x0=np.zeros(cir.n - 1), maxiterations=100,
                          x0_unknown=False)
            assert pss.converged
            K0 = PAC(cir).covariance(pss)
            ib = [str(nd) for nd in cir.nodes].index('b')
            ratios.append(float(K0[ib, ib]) / (kT / Cc))
        assert ratios[-1] > 0.95 and ratios[-1] < 1.05, \
            '%s: covariance is %.4f of kT/C at 1600 points' % (method, ratios[-1])
        assert abs(ratios[-1] - 1.0) < abs(ratios[0] - 1.0), \
            '%s: not converging toward kT/C (%s)' % (method, ratios)


def test_a_coloured_covariance_takes_a_modulated_non_power_law_source():
    """Coloured noise that is not a power law AND follows the orbit
    (2026-09-26; Andreas: "Do 2a and 2b").  SEPARABLE -- a level that follows
    the state under a fixed spectral shape, ``C(x, w) = C(x, w_ref) s(w)``
    -- replays one amplitude per point and weights each band frequency by
    `s`; with the SHAPE moving along the orbit the density is read at every
    point for every band frequency (the quasi-static model `pnoise` and
    `sampled_variance` use), warned for its cost and its sign-blind root.

    Measured (radau, a driven RC, the level from a positive clock):
      * separable, against the same noise realised as a white source
        through an explicit RC filter times the clock (`_NuMult`) -- the
        white Lyapunov route: -3.6e-6 at t = 0, 3.6e-6 at worst over a
        profile that varies 2.8x in the period (the band below `fmin`);
      * the separable source forced through the moving-shape path: 2.2e-16;
      * the corner moving (`shape = 0.5`) against `sampled_variance` (the
        adjoint route, the same quasi-static model): +6.0e-5 / +6.4e-5 at
        two instants -- its own quadrature and the holes its fold leaves
        (not asserted: that reference costs 235 s);
      * the moving-shape path at ``shape -> 0`` meets the separable one."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    Rf, Cf, g, Pw = 1e3, 0.3e-9, 1e-3, 1e-20
    tau, P = Rf * Cf, g * g * Pw * Rf * Rf

    def build(kind, shape=0.0, npts=200):
        c = SubCircuit()
        for nd in ('lo', 'out'):
            c.add_node(nd)
        c['Vlo'] = VSin('lo', gnd, va=1.0, vo=1.5, freq=1.0 / T)
        c['Ro'] = R('out', gnd, r=1e3, noisy=False)
        c['Co'] = C('out', gnd, c=0.5e-9)
        if kind == 'element':
            c['n'] = _ModLorentzCtl('out', gnd, 'lo', gnd, noisePSD=P, tau=tau,
                                    k=1.0, shape=shape)
        else:
            c.add_node('f')
            c['nw'] = IS('f', gnd, i=0.0, noisePSD=Pw)
            c['rf'] = R('f', gnd, r=Rf, noisy=False)
            c['cf'] = C('f', gnd, c=Cf)
            c['mx'] = _NuMult('out', gnd, 'f', gnd, 'lo', gnd, k=g)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, maxiterations=40)
        o = [str(x) for x in c.nodes if str(x) != 'gnd!'].index('out')
        return pss, o, PAC(c, toolkit=circuit.numeric)

    fmin = 1e-6 / T
    ## separable, against its white-through-filter realisation
    pss, o, pac = build('element')
    with pytest.warns(RuntimeWarning, match='SIGN-BLIND'):
        _K, se = pac.covariance(pss, samples=True, colour_fmin=fmin)
    pf, of, pacf = build('filtered')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _Kf, sf = pacf.covariance(pf, samples=True)
    a = np.array([K[o, o] for K in se])
    b = np.array([K[of, of] for K in sf])
    n = min(len(a), len(b))
    assert np.max(np.abs(a[:n] / b[:n] - 1.0)) < 1e-5, np.max(np.abs(a[:n] / b[:n] - 1.0))
    ## the classifier, and the two paths on one source
    host = pss._lyapunov_host()
    fp = host._state_map()
    _cnt, states = pac._injection_points(host, fp)
    wsp = [2 * np.pi * f for f in (1e3, 1e6, 1e7, 5e7)]
    assert separable([pac._noise_components(host, states).one_element_cy(('n',), w)
                      for w in wsp])
    p5, o5, pac5 = build('element', shape=0.5, npts=100)
    h5 = p5._lyapunov_host()
    _c5, st5 = pac5._injection_points(h5, h5._state_map())
    assert not separable([pac5._noise_components(h5, st5).one_element_cy(('n',), w)
                          for w in wsp])
    p0, o0, pac0 = build('element', npts=100)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Ks = pac0.covariance(p0, colour_fmin=fmin)[o0, o0]
        ## (on the instance's factory: the class keeps its own)
        consulted = []
        _noise_seam(pac0, separable=staticmethod(
            lambda Cs, tol=1e-9: consulted.append(1) or False))
        Kn = pac0.covariance(p0, colour_fmin=fmin)[o0, o0]
    ## ⚠ the equality below holds VACUOUSLY if the patch is never looked up
    ## (a refactor that stops reaching `separable` through the factory)
    assert consulted, 'the patched separable was never consulted'
    assert abs(Kn / Ks - 1.0) < 1e-12, Kn / Ks - 1.0
    ## the moving-shape path meets the separable one as the corner stops
    pt, ot, pact = build('element', shape=1e-6, npts=100)
    with pytest.warns(RuntimeWarning, match='SHAPE changes along the orbit'):
        Kt = pact.covariance(pt, colour_fmin=fmin)[ot, ot]
    assert abs(Kt / Ks - 1.0) < 1e-4, Kt / Ks - 1.0


def test_the_coloured_band_integral_resolves_a_high_q_line():
    """The band integral's grid is ADAPTIVE (2026-09-25).  A driven parallel
    tank at Q = 20, resonant at 1.37 f0, with a 1/f current: on the fixed
    40-per-decade log grid the coloured covariance read -10.5 % against 640
    per decade (-0.16 % at 160; Q = 5: -1.4e-4) -- the line is 2.5 % wide
    and the grid 6 %.  Adaptive Simpson on the log axis (a panel against
    itself on five points, ``|S_2 - S_1| / 15``, on its share of
    `COLOURED_REFINE_TOL`): 10 and 40 per decade agree to 1.3e-7."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    fr = 1.37 / T
    Lv = 1e-6
    Cv = 1.0 / ((2.0 * np.pi * fr) ** 2 * Lv)
    Rp = 20.0 * 2.0 * np.pi * fr * Lv
    c = SubCircuit()
    c.add_node('in')
    c.add_node('out')
    c['V'] = VSin('in', gnd, va=0.1, vo=0.0, freq=1.0 / T)
    c['Rs'] = R('in', 'out', r=10.0 * Rp, noisy=False)
    c['Rp'] = R('out', gnd, r=Rp, noisy=False)
    c['C'] = C('out', gnd, c=Cv)
    c['L'] = L('out', gnd, L=Lv)
    c['n'] = _Flicker('out', gnd, i=0.0, noisePSD=1e-20, fref=1.0)
    pss = PSS(c, method='radau', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 100, maxiterations=40)
    o = [str(x) for x in c.nodes if str(x) != 'gnd!'].index('out')
    pac = PAC(c, toolkit=circuit.numeric)
    a = pac.covariance(pss, colour_fmin=1e3, points_per_decade=10)[o, o]
    b = pac.covariance(pss, colour_fmin=1e3, points_per_decade=40)[o, o]
    assert abs(a / b - 1.0) < 1e-6, a / b - 1.0


def test_a_switched_capacitor_holds_kTC_with_per_step_CY():
    """✅ AN EXTERNAL CROSS-CHECK'S FINDING, CLOSED: `covariance` evaluated one
    `CY` for the whole period, so a switch's `4kT g(t)` was outside its
    FORMULATION (the accumulation was already per step), and the
    cyclostationarity refusal in `_cy_reduced` closed the door on the
    smallest circuit with `kT/C` noise.  With `CY` per step at the step's
    own state, the sample-and-hold (`Ron = 1k`, `Roff = 1G`, `C = 100 pF`,
    `f = 100 kHz`, the suite's fixture to the parameter) gives, over the
    period:

        points     200      400      800      1600
        hold       0.9920   0.9987   0.9998   1.0000     x kT/C   (second order)
        track      0.7398   0.8467   0.9157   0.9556     x kT/C   (the known O(h/tau))

    A reference simulator's sampled pnoise reads 0.99999 x kT/C at every instant of the
    hold phase.  The tracking phase sits at the covariance's O(h/tau)
    floor, IDENTICAL to the switch-held-closed control (no clock swing) at
    every grid -- that floor is a separate, recorded item and not this
    one.  Gated on the hold-phase mean and the control identity at 400
    points, and on the hold phase's order.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    fclk, cval, kb, temp = 100e3, 100e-12, 1.38e-23, 300.0
    ktc = kb * temp / cval
    T = 1.0 / fclk

    def build(vth, vck):
        cir = SubCircuit()
        cir.add_node('in')
        cir.add_node('out')
        cir.add_node('ck')
        cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
        cir['Vck'] = VSin('ck', gnd, vo=0.0, va=vck, freq=fclk, phase=90.0)
        cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                               vth=vth, vs=50e-3, temp=temp, kb=kb)
        cir['C0'] = C('out', gnd, c=cval)
        return cir

    def run(vth, vck, npts):
        cir = build(vth, vck)
        pss = PSS(cir, method='gear', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                      maxiterations=100)
        assert pss.converged
        io = [str(n) for n in cir.nodes if str(n) != 'gnd!'].index('out')
        K0, Ks = PAC(cir, toolkit=circuit.numeric).covariance(pss,
                                                              samples=True)
        v = np.array([np.asarray(k, dtype=float)[io, io] for k in Ks]) / ktc
        tt = np.linspace(0.0, T, len(v), endpoint=False)
        hold = v[(tt > 0.3 * T) & (tt < 0.45 * T)].mean()
        track = v[tt < 0.2 * T].mean()
        return float(np.asarray(K0, dtype=float)[io, io]) / ktc, hold, track
    k0_sw, hold4, track4 = run(0.0, 1.0, 400)
    k0_ctrl, _h, track_ctrl = run(-1.0, 0.0, 400)
    assert abs(hold4 - 1.0) < 3e-3, \
        'held variance %.4f x kT/C at 400 points; the reference samples 0.99999' \
        % hold4
    assert abs(track4 / track_ctrl - 1.0) < 1e-6, \
        'the tracking phase (%.4f) must sit on the control (%.4f): the same ' \
        'O(h/tau) floor' % (track4, track_ctrl)
    assert abs(k0_sw / k0_ctrl - 1.0) < 1e-6
    _k, hold8, _t = run(0.0, 1.0, 800)
    assert abs(hold8 - 1.0) < 8e-4 and (1.0 - hold4) / (1.0 - hold8) > 3.0, \
        'the held variance must converge at second order: %.2e -> %.2e' \
        % (1.0 - hold4, 1.0 - hold8)


def test_trbdf2_covariance_converges_to_kTC_at_second_order():
    """TR-BDF2's OWN per-step injection (DAE-projected Van Loan) makes the
    covariance converge to the exact `kT/C` at SECOND order -- against
    Gear-2's first.

    `Var(v_C) = kT/C` is exact and independent of R (see
    `test_the_periodic_covariance_converges_to_kTC`).  Gear-2's covariance
    reaches it at O(h) (a piecewise-constant white-noise injection); the
    TR-BDF2 Van Loan injection reaches it at O(h^2), so at a coarse grid it
    is more than two orders closer -- 3.8e-4 vs 6.9e-2 at 100 points here.
    This is the injection surface's payoff, the same direction as the
    monodromy/PPV gains.

    ⚠ THE ASSERTION IS THE RATE (~4x per doubling), NOT machine precision.
    A machine-zero kT/C would mean a method-consistent `Q = P(1 - A^2)`
    fudge that makes the stationary observable exact by absorbing the
    method's error -- and gives the WRONG transient covariance.  The error
    must be O(h^2) and DECREASING, which is what pins the injection as the
    exact per-step integral rather than a stationary fit.
    """
    import warnings
    from pycircuit.circuit.constants import kboltzmann
    circuit.default_toolkit = circuit.numeric

    def ratio(npts, method, Cval=1e-7, per=1e-3):
        cir = _rc_noisy(Cval=Cval, per=per)
        pss = PSS(cir, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=per, timestep=per / npts, maxiterations=40)
        assert pss.converged
        K0 = PAC(cir).covariance(pss)
        irn = pss.irefnode
        k = cir.get_node_index(cir.get_node('b'))
        k = k - 1 if k > irn else k
        T = float(circuit.defaultepar.T)
        return K0[k, k] / (kboltzmann * T / Cval)

    errs = [abs(1.0 - ratio(n, 'trbdf2')) for n in (100, 200, 400)]
    ## second order: ~4x per doubling (allow a band)
    for a, b in zip(errs, errs[1:]):
        assert 3.2 < a / b < 5.0, \
            'trbdf2 covariance falls %.2fx per doubling (%s), not the ~4x ' \
            'of O(h^2)' % (a / b, errs)
    ## and it is already tight at the coarsest grid, far better than gear's
    assert errs[0] < 2e-3, errs
    gear0 = abs(1.0 - ratio(100, 'gear'))
    assert gear0 / errs[0] > 20.0, \
        'trbdf2 covariance should be >20x closer than gear at 100 pts; ' \
        'gear %.2e vs trbdf2 %.2e' % (gear0, errs[0])


def test_radau_covariance_converges_to_kTC_faster_than_second_order():
    """Radau IIA(3)'s covariance converges to the exact `kT/C` FASTER than
    second order -- the payoff of pairing an exact injection with an order-5
    transition.

    The per-step injection `Q_n` is the SAME DAE-projected Van Loan integral
    TR-BDF2 uses -- the exact continuous per-step covariance, method-
    independent.  On this LTI RC fixture the operating point is constant, so
    that injection is exact and the ONLY discretisation error left in the
    Lyapunov recursion `K = A K A^T + Q` is the transition `A_n`'s departure
    from `exp(A h)`: O(h^2) for TR-BDF2, O(h^5) for Radau.  So Radau's
    covariance falls far faster than TR-BDF2's ~4x per doubling -- measured
    ~30x -- and reaches `kT/C` to ~1e-9 at 100 points.

    ⚠ THE RATE IS SET BY THE TRANSITION, NOT THE INJECTION, ON AN LTI ORBIT.
    On a time-varying orbit the injection's single-point linearisation per
    step would become the bottleneck and the rate would drop back toward the
    injection's order; this fixture isolates the transition, which is the
    point of the comparison with TR-BDF2 on the identical fixture.  The
    assertion is still the RATE and the closeness, never machine zero (a
    stationary fit that hit kT/C exactly would corrupt the transient
    covariance -- see `test_trbdf2_covariance_converges_to_kTC`).
    """
    import warnings
    from pycircuit.circuit.constants import kboltzmann
    circuit.default_toolkit = circuit.numeric

    def ratio(npts, method, Cval=1e-7, per=1e-3):
        cir = _rc_noisy(Cval=Cval, per=per)
        pss = PSS(cir, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=per, timestep=per / npts, maxiterations=40)
        assert pss.converged
        K0 = PAC(cir).covariance(pss)
        irn = pss.irefnode
        k = cir.get_node_index(cir.get_node('b'))
        k = k - 1 if k > irn else k
        T = float(circuit.defaultepar.T)
        return K0[k, k] / (kboltzmann * T / Cval)

    errs = [abs(1.0 - ratio(n, 'radau')) for n in (100, 200, 400)]
    ## faster than second order: each doubling cuts the error by far more than
    ## the ~4x of O(h^2) (the exact injection lets the order-5 transition show)
    for a, b in zip(errs, errs[1:]):
        assert a / b > 8.0, \
            'radau covariance falls only %.2fx per doubling (%s); the exact ' \
            'injection + order-5 transition should beat O(h^2)' % (a / b, errs)
    ## and it is already at the monodromy floor at the coarsest grid
    assert errs[0] < 1e-6, errs
    ## far closer than gear's first-order injection at the same grid
    gear0 = abs(1.0 - ratio(100, 'gear'))
    assert gear0 / errs[0] > 20.0, \
        'radau covariance should be >20x closer than gear at 100 pts; ' \
        'gear %.2e vs radau %.2e' % (gear0, errs[0])


def test_the_stage_method_covariance_injects_at_the_stages_and_holds_kTC_across_a_switching_edge():
    """⚠⚠ `covariance` under radau/trbdf2 read a switch's HELD variance O(h)
    low (2026-09-15): its Van Loan injection froze `C`, `G`, `CY` at the END of
    each step, i.e. at the OFF conductance across the switch-off edge --
    1 - 0.876 / 0.934 / 0.966 kT/C at 400 / 800 / 1600 points under radau,
    0.87 under trbdf2 -- while gear read 0.9987.  Now the source enters every
    stage (`_stage_injection`): measured radau 1.3e-7 off held and tracking at
    400 points, trbdf2 2.5e-3 held at 400 and 1.0e-2 at 200 (second order).
    The constant-operating-point convergence tests beside this one keep their
    rates (radau 1.4e-9 at 100 points, trbdf2 ~4x per doubling)."""
    ktc = _KB * _TEMP / 100e-12
    err = {}
    for method, npts in (('radau', 400), ('trbdf2', 200), ('trbdf2', 400)):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), method, npts)
        N = len(pss.factored_period().steps)
        _K0, Ks = pac.covariance(pss, samples=True)
        held = float(np.asarray(Ks[int(0.375 * N)], dtype=float)[io, io]) / ktc
        track = float(np.asarray(Ks[int(0.1 * N)], dtype=float)[io, io]) / ktc
        err[(method, npts)] = (1.0 - held, 1.0 - track)
    assert abs(err[('radau', 400)][0]) < 1e-6, err
    assert abs(err[('radau', 400)][1]) < 1e-6, err
    assert abs(err[('trbdf2', 400)][0]) < 5e-3, err
    r = err[('trbdf2', 200)][0] / err[('trbdf2', 400)][0]
    assert 3.2 < r < 5.0, (r, err)


def test_the_esdirk43_covariance_holds_kTC_across_a_switching_edge_despite_its_negative_weight():
    """ESDIRK43's weights are `0.158, 0, 0.187, 0.681, -0.275, 0.25`, so the
    per-stage injection radau and trbdf2 use (variance `CY/(2 h b_i)`) is
    undefined, and it kept the END-of-step Van Loan: held variance
    1 - 0.367 / 0.224 / 0.124 kT/C at 100 / 200 / 400 points on the sampler.
    It now takes the Van Loan injection at the stage states with positive
    trapezoid weights over the abscissae: measured 2.1e-2 / 3.4e-3 / 3.7e-4,
    tracking unchanged (2.6e-6 at 400), and identical to the end-of-step form
    wherever the operating point is constant (the weights sum to one)."""
    ktc = _KB * _TEMP / 100e-12
    err = {}
    for npts in (200, 400):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), 'esdirk43', npts)
        N = len(pss.factored_period().steps)
        _K0, Ks = pac.covariance(pss, samples=True)
        held = float(np.asarray(Ks[int(0.375 * N)], dtype=float)[io, io]) / ktc
        track = float(np.asarray(Ks[int(0.1 * N)], dtype=float)[io, io]) / ktc
        err[npts] = (1.0 - held, 1.0 - track)
    assert abs(err[400][0]) < 1e-3, err
    assert err[200][0] / err[400][0] > 4.0, err
    assert abs(err[400][1]) < 1e-4, err


def test_the_driven_fold_is_measured_for_accuracy_on_radau_not_only_alignment():
    """Item 5 of the non-uniform-grid list (2026-09-21): the fold on a
    DRIVEN circuit, measured.  Sine-clocked sampler, radau, against a
    radau uniform-3200 reference: the 95-point fold reads 3.7e-5 of the
    output swing where uniform 95 reads 2.4e-4 and uniform 190 5.1e-5;
    held variance 1.00000 x kT/C on the fold, 0.99992 / 0.99999 uniform.
    (trbdf2 261: 1.5e-4 vs 2.1e-3; gear 125: 1.2e-3 vs 4.5e-2 -- but
    gear's held variance on its fold is 0.78 against 0.98 uniform: the
    fold resolves the STATE, and gear's covariance has its recorded
    O(h/tau) tracking floor where the state is flat -- see `lte_grid`.)
    Pinned: the fold at least 4x better than uniform at its own count on
    the waveform, held within 5e-4 of kT/C.

    ⚠ 2026-09-28 the stage methods honour `relref`, and `lte_grid`'s
    pre-run takes its own, default 'pointlocal': on 'sigglobal' this fold
    had 54 points and held 0.99950 x kT/C, no better than uniform (the
    switch opening left coarse); on 'pointlocal' 92 points, 3.7e-5 of the
    swing, held 0.999996.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    fclk, cval, kb, temp = 100e3, 100e-12, 1.38e-23, 300.0
    ktc = kb * temp / cval
    T = 1.0 / fclk

    def build():
        cir = SubCircuit()
        cir.add_node('in')
        cir.add_node('out')
        cir.add_node('ck')
        cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0)
        cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
        cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                               vth=0.0, vs=50e-3, temp=temp, kb=kb)
        cir['C0'] = C('out', gnd, c=cval)
        return cir

    def out_of(cir, p):
        io = [str(n_) for n_ in cir.nodes].index('out')
        io = io if io < p.irefnode else io - 1
        X = np.asarray(p.waveform[1], float)
        X = X if X.shape[0] == cir.n - 1 else np.delete(X, p.irefnode, axis=0)
        return io, np.asarray(p.waveform[0], float), X[io]

    def solve(npts, grid=None, seed=None, reltol=1e-8):
        cir = build()
        p = PSS(cir, method='radau', reltol=reltol)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / npts, grid=grid,
                    x0=np.zeros(cir.n - 1) if seed is None else seed,
                    maxiterations=100)
        assert p.converged
        return cir, p

    ## ⚠ 3200, not 1600: a 1600-point radau reference has ~1e-4 of the
    ## swing left in it and reads the fold's 3.7e-5 as 1.1e-4 (ratio 2.4)
    cr, pr = solve(3200, reltol=1e-10)
    _io, tsr, vr = out_of(cr, pr)
    swing = vr.max() - vr.min()
    cir = build()
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        fr, seed = PSS(cir, method='radau', reltol=1e-8).lte_grid(
            T, x0=np.zeros(cir.n), reltol=1e-5)
    fr = np.asarray(fr, float)
    errs, helds = {}, {}
    for label, grid, sd in (('fold', fr, seed), ('uniform', None, None)):
        c2, p2 = solve(len(fr), grid=grid, seed=sd)
        io, ts, v = out_of(c2, p2)
        errs[label] = float(np.max(np.abs(v - np.interp(ts % T, tsr, vr))) / swing)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _K0, Ks = PAC(c2, toolkit=circuit.numeric).covariance(p2, samples=True)
        vv = np.array([np.asarray(k, float)[io, io] for k in Ks]) / ktc
        tt = ts[:len(vv)]
        helds[label] = float(vv[(tt > 0.3 * T) & (tt < 0.45 * T)].mean())
    assert errs['uniform'] / errs['fold'] > 4.0, errs
    assert abs(helds['fold'] - 1.0) < 5e-4, helds


def test_covariance_on_a_staged_solve_borders_its_lyapunov_closure_with_the_moving_events():
    """Events phase B (2026-09-22): `PAC.covariance` on a staged solve
    closes on the TOTAL monodromy with the noise-driven motion of the
    crossings in the injection, and samples at FIXED times.

    Measured on `_jitter_sampler` at 100 / 200 / 400 points (radau, the
    tanh switch): bordered held variance 0.9991 / 0.9992 / 0.9992 of the
    analytic `(s2/s1)^2 kT/C_n` (flat in the count -- the residual is the
    100 ps tracking lag of the fixture, not the grid); the UNBORDERED
    closure on the same staged solve read 4.08 / 3.79 / 3.29.  ⚠ With
    VSwitch's COMPACT transition (2026-09-22) both read 0.9996 at 100
    points: the 4x was the tanh's tails outside the landed window, which
    the per-step maps could not resolve; with the transition inside the
    window they carry the threshold's motion themselves and the bordering
    corrects a sliver.  And the node-rate
    correction is not decoration: without it the NOISELESS sawtooth
    source reads 0.225 kT/C_n at t = 0.7T -- exactly its node's share of
    the moving segment, ((0.925 - 0.7) / (0.925 - 0.45))^2."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    cir = _jitter_sampler(T)
    pss = PSS(cir, method='radau', reltol=1e-9)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 100, maxiterations=100,
                  state_events=True)
    assert pss.converged
    ## the turn-off crossing lands at 0.45 T (saw = 2.5 V), both edges
    fr = sorted(float(e) for e in pss.event_times if 0.1 < e < 0.9)
    assert len(fr) == 2 and abs(fr[0] - 0.45) < 1e-4 and abs(fr[1] - 0.45) < 1e-4
    names = [str(n) for n in cir.nodes]
    red = [n for i, n in enumerate(names) if i != pss.irefnode]
    ih, inn, isaw = red.index('hold'), red.index('n'), red.index('saw')
    ts = np.asarray(pss.waveform[0], dtype=float)
    kT = 1.380649e-23 * 300.0
    exp_n = kT / 1e-12
    exp_h = 4.0 * exp_n
    pac = PAC(cir, toolkit=circuit.numeric)
    K0, seq = pac.covariance(pss, samples=True)
    assert len(seq) == len(ts)
    assert np.linalg.norm(seq[-1] - K0) / np.linalg.norm(K0) < 1e-8
    var_n = np.mean([K[inn, inn] for K in seq])
    assert abs(var_n / exp_n - 1.0) < 5e-3, var_n / exp_n
    held = seq[int(np.searchsorted(ts, 0.7 * T))][ih, ih]
    assert abs(held / exp_h - 1.0) < 5e-3, held / exp_h
    ## the hold is flat while held, silent while tracking
    for ph in (0.5, 0.6, 0.8):
        assert abs(seq[int(np.searchsorted(ts, ph * T))][ih, ih] / held - 1.0) < 1e-3
    assert seq[int(np.searchsorted(ts, 0.2 * T))][ih, ih] < 1e-6 * exp_h
    ## a noiseless source carries no variance -- the samples are at fixed times
    j7 = int(np.searchsorted(ts, 0.7 * T))
    assert seq[j7][isaw, isaw] < 1e-9 * exp_n, seq[j7][isaw, isaw] / exp_n
    ## the unbordered closure on the same solve is wrong by O(1)
    ev = pss._event_columns
    pss._event_columns = None
    try:
        _K0u, sequ = pac.covariance(pss, samples=True)
    finally:
        pss._event_columns = ev
    ## ⚠ RE-PINNED 2026-09-22 for VSwitch's COMPACT transition: the unbordered
    ## closure reads 0.9996 too -- with the whole transition inside the
    ## landed window the per-step maps carry the threshold's motion.  The
    ## 4.08 / 3.79 / 3.29 it read with the tanh was the tails.
    assert abs(sequ[j7][ih, ih] / exp_h - 1.0) < 5e-3, sequ[j7][ih, ih] / exp_h
    ## and the node-rate correction is what keeps the source silent
    pac._orbit_rate = lambda p, nodes: np.zeros((len(p.waveform[0]), p.cir.n - 1))
    try:
        _K0z, seqz = pac.covariance(pss, samples=True)
    finally:
        del pac._orbit_rate
    share = ((0.925 - ts[j7] / T) / (0.925 - 0.45)) ** 2
    assert abs(seqz[j7][isaw, isaw] / exp_n / share - 1.0) < 2e-2, (seqz[j7][isaw, isaw] / exp_n, share)


def test_event_jitter_is_the_crossings_own_noise_and_matches_the_analytic_sigma():
    """The crossings' noise-driven jitter as a user-facing quantity
    (2026-09-23): `PAC.event_jitter` reads it off the bordered Lyapunov
    closure -- ``Cov(dtheta) = dth K_0 dth^T + Gt^-1 D Gt^-T``, the
    stationary state at the period start plus this period's per-step
    injections.  On `_jitter_sampler` (a sawtooth of slope `s_1 = V1/0.9T`
    crossing a threshold node carrying `kT/C_n`) the analytic answer is
    ``sqrt(kT/C_n) / s_1`` at the turn-off crossing: 11.5844 ps.  Measured
    11.5873 ps -- 1.0002 -- and FLAT at 100 / 200 / 400 points, the same
    crossing motion that gives the held capacitor its `(s_2/s_1)^2 kT/C_n`
    in `covariance`.  The reset edges read 0.6437 ps (the threshold's own
    faster slope there).  Pinned: within 1 % of the analytic sigma at 100
    and 200 points, the two edges of one window equal, and an oscillator
    refused."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    kT = 1.380649e-23 * 300.0
    sigma_exact = np.sqrt(kT / 1e-12) / (5.0 / (0.9 * T))
    out = {}
    for N in (100, 200):
        cir = _jitter_sampler(T)
        pss = PSS(cir, method='radau', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / N, maxiterations=100, state_events=True)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            j = pac.event_jitter(pss)
        assert len(j['sigma_t']) == 4 and len(j['fractions']) == 4
        fr = np.asarray(j['fractions'], dtype=float)
        assert abs(fr[0] - 0.45) < 1e-4 and abs(fr[2] - 0.925) < 1e-3, fr
        ## the two edges of one window cross together: the same sigma
        assert abs(j['sigma_t'][0] / j['sigma_t'][1] - 1.0) < 1e-6
        assert abs(j['sigma_t'][2] / j['sigma_t'][3] - 1.0) < 1e-6
        ## the reset edge is faster, so it jitters less
        assert j['sigma_t'][2] < 0.1 * j['sigma_t'][0]
        out[N] = float(j['sigma_t'][0])
    for N, s in out.items():
        assert abs(s / sigma_exact - 1.0) < 1e-2, (N, s, sigma_exact)
    assert abs(out[100] / out[200] - 1.0) < 1e-3, out
    ## an oscillator's crossings diffuse with its phase: refused, not fudged
    osc = _comparator_relaxation_oscillator()
    seed, To = _relaxation_oscillator_seed(osc)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        qo = PSS(osc, method='radau', reltol=1e-9)
        qo.solve(period=To, timestep=To / 200,
                 x0=seed, maxiterations=100, state_events=True)
    try:
        PAC(osc, toolkit=circuit.numeric).event_jitter(qo)
        assert False, 'an oscillator must be refused'
    except ValueError as e:
        assert 'diffuse' in str(e)


def _rc_flicker(method, npts, T=1e-6, Rv=1e3, Cv=1e-9, k=1e-20):
    """A driven LTI RC (`fc = 159 kHz`, `f0 = 1 MHz`) whose ONLY noise is a
    1/f current `k/f` into the capacitor node -- the closed-form gate of a
    coloured covariance."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('in')
    c.add_node('out')
    c['V'] = VSin('in', gnd, va=0.1, vo=0.0, freq=1.0 / T)
    c['R'] = R('in', 'out', r=Rv, noisy=False)
    c['C'] = C('out', gnd, c=Cv)
    c['n'] = _Flicker('out', gnd, i=0.0, noisePSD=k, fref=1.0)
    pss = PSS(c, method=method, reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, maxiterations=40)
    o = [str(x) for x in c.nodes if str(x) != 'gnd!'].index('out')
    return c, pss, o, PAC(c, toolkit=circuit.numeric)


def _rc_band_variance(f1, f2, Rv=1e3, Cv=1e-9, k=1e-20):
    """``int_f1^f2 (k/f) R^2 / (1 + (f/fc)^2) df`` in closed form."""
    fc = 1.0 / (2.0 * np.pi * Rv * Cv)
    F = lambda f: 0.5 * np.log(f * f / (f * f + fc * fc))     # noqa: E731
    return k * Rv ** 2 * (F(f2) - F(f1))


def test_a_coloured_covariance_integrates_the_band_against_the_closed_form():
    """Coloured noise in `covariance` (2026-09-25): the white part of every
    source goes through the Lyapunov recursion, the coloured part is
    integrated in the FREQUENCY DOMAIN over the band `[fmin, fmax]` --
    per input frequency the forced periodic response to the MODULATED
    source, ``K = int (w1/w)^EF Re[y y^H] dnu`` -- no shaping filter, no
    fitted Lorentzians (`PAC._coloured_covariance`).

    Gate: an LTI RC with only a 1/f current, against ``int (k/f) |H|^2 df``
    in closed form, over `[1 kHz, the grid's Nyquist]`.  Measured (rel.):

        N      gear       trbdf2     radau / glm2 (its twin)   euler
        100   -5.5e-5    -8.5e-6    -5.5e-11                   -9.8e-4
        200   -1.6e-5    -2.4e-6    -3.3e-11                   -4.9e-4
        400   -4.5e-6    -6.6e-7    -8.2e-12                   -2.5e-4

    -- each method's discrete transfer at its own order.  The band
    quadrature (`_power_law_weights`: the power law exact per interval,
    Richardson on top) is fourth order: 5.9e-6 / 8.3e-10 / 2.0e-11 at 5 /
    10 / 20 per decade.  (Until the Richardson step, 2026-09-25, radau read
    a flat +4e-9 here -- the trapezoid in ln nu, second order; in `nu` it
    would carry 6e-4 of a 1/f band at 40 per decade.)  The samples of a
    stationary answer are flat (1e-15), the first is `K0`.
    Refused: no `fmin` (a 1/f variance grows as ln(fmax/fmin) without
    limit), a band past the grid's Nyquist, and the plain trapezoidal
    map, whose covariance lives on the pair (x, iq)."""
    fmin = 1e3
    for method, N, tol in (('radau', 100, 1e-9), ('gear', 100, 1e-4),
                           ('gear', 200, 3e-5)):
        _c, pss, o, pac = _rc_flicker(method, N)
        Ns = len(pss.factored_period().steps)
        fmax = 0.5 * Ns / float(pss.period)
        K0, seq = pac.covariance(pss, samples=True, colour_fmin=fmin)
        rel = K0[o, o] / _rc_band_variance(fmin, fmax) - 1.0
        assert abs(rel) < tol, (method, N, rel)
        assert np.array_equal(K0, seq[0]) and len(seq) == Ns + 1
        v = np.array([Kj[o, o] for Kj in seq])
        assert (v.max() - v.min()) < 1e-12 * v.mean(), (method, N)
        if method == 'gear' and N == 100:
            rel100 = rel
    ## second order: 5.5e-5 -> 1.6e-5
    assert 2.5 < rel100 / rel < 5.0, (rel100, rel)
    ## an inner band, fmax below the Nyquist
    _c, pss, o, pac = _rc_flicker('radau', 100)
    K = pac.covariance(pss, colour_fmin=1e2, colour_fmax=1e7)
    assert abs(K[o, o] / _rc_band_variance(1e2, 1e7) - 1.0) < 1e-7
    ## (a TypeError, a required argument missing, since 2026-09-29)
    with pytest.raises(TypeError,
                       match='ln\\(colour_fmax/colour_fmin\\)'):
        pac.covariance(pss)
    with pytest.raises(ValueError, match='Nyquist'):
        pac.covariance(pss, colour_fmin=1e3, colour_fmax=1e9)
    _c, pss, o, pac = _rc_flicker('trap', 100)
    with pytest.raises(NotImplementedError, match='pair \\(x, iq\\)'):
        pac.covariance(pss, colour_fmin=fmin)


def test_event_jitter_integrates_a_coloured_threshold_and_keeps_the_white_part():
    """Coloured noise in `event_jitter` (2026-09-25): `_jitter_sampler` plus a
    1/f current `k/f` into the threshold node `n` (`R_n C_n`, LTI).  The
    crossing moves by ``-v_n/s_1``, so its coloured variance is exactly
    ``Var_band(v_n) / s_1^2``, `Var_band` the RC's closed form over `[fmin,
    fmax]`; the held capacitor takes `(s_2/s_1)^2` of it.  Measured at 100
    and 200 points (radau): the coloured part of sigma^2 0.999999 of
    exact, the node's coloured variance 0.999999, the hold's 0.99995 (the
    fixture's tracking lag, as for the white part); gear and trbdf2
    0.99999 on the jitter.

    ⚠ THE WHITE PART MUST NOT SEE THE FLICKER.  The Lyapunov pieces read
    `CY` at `w0`, which on this circuit carries the flicker at 1 MHz; they
    read the component model's WHITE part while a coloured covariance runs
    (`PAC._lyap_cy`).  The white-only run is the reference: the flicker
    run minus the closed-form coloured part returns it to 1e-6."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    Rn, Cn, k = 1e4, 1e-12, 4e-18
    fmin, fmax = 1e2, 2e6
    s1 = 5.0 / (0.9 * T)
    var_n = _rc_band_variance(fmin, fmax, Rv=Rn, Cv=Cn, k=k)
    out = {}
    for flick in (False, True):
        cir = _jitter_sampler(T)
        if flick:
            cir['Fn'] = _Flicker('n', gnd, i=0.0, noisePSD=k, fref=1.0)
        pss = PSS(cir, method='radau', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=100,
                      state_events=True)
        assert pss.converged
        pac = PAC(cir, toolkit=circuit.numeric)
        band = dict(colour_fmin=fmin, colour_fmax=fmax) if flick else {}
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            j = pac.event_jitter(pss, **band)
            K0, seq = pac.covariance(pss, samples=True, **band)
        red = [str(x) for i, x in enumerate(cir.nodes) if i != pss.irefnode]
        tms = np.asarray(pss.factored_period().times, dtype=float)
        jh = int(np.argmin(abs(tms - 0.7 * T)))
        out[flick] = (float(j['sigma_t'][0]) ** 2,
                      seq[jh][red.index('n'), red.index('n')],
                      seq[jh][red.index('hold'), red.index('hold')])
        if flick:
            with pytest.raises(TypeError, match='COLOURED'):
                pac.event_jitter(pss)
    (sw, nw, hw), (sc, nc, hc) = out[False], out[True]
    assert abs((sc - sw) / (var_n / s1 ** 2) - 1.0) < 1e-5
    assert abs((nc - nw) / var_n - 1.0) < 1e-5
    assert abs((hc - hw) / (4.0 * var_n) - 1.0) < 2e-4


def test_a_coloured_covariance_takes_a_flicker_whose_exponent_differs_between_entries():
    """A coloured component whose power-law exponent differs between its
    entries was refused by the band integral ("no one amplitude to
    replay") until 2026-09-26.  Now:

      * entries in DISJOINT index blocks, one exponent each (sources of
        different slope on branches sharing no node) split EXACTLY into
        independent components (`_split_by_exponent`): identical to the
        same sources as two elements, and to the closed form -3.2e-10 /
        +7.4e-11 (1/f^0.8 / 1/f^2);
      * correlated entries of different slope cannot split: the density
        `B (w1/w)^EF` is rooted at every point per band frequency (the
        moving-shape path), warned; against the closed form -1.0e-8 /
        -4.6e-8 / -2.2e-8 (both diagonals and the cross entry).
    (Two sources SHARING a node sum in one entry, which the component
    model cannot fit as one power law: that element is per-band already.)
    """
    from scipy.integrate import quad
    import warnings as _w
    K = {}
    for kind in ('mixed', 'split', 'corr'):
        pss, pac, names, (T, Rv, Cv, k) = _mixed_exponent_rc(kind)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            K[kind] = pac.covariance(pss, colour_fmin=1e3)
        if kind == 'corr':
            assert any('different power-law exponents' in str(r.message)
                       for r in rec)
    fmax = 0.5 * len(pss.factored_period().steps) / T
    fc = 1.0 / (2.0 * np.pi * Rv * Cv)

    def ref(ef):
        return quad(lambda u: k * np.exp((1.0 - ef) * u) * Rv ** 2
                    / (1.0 + np.exp(2.0 * u) / fc ** 2),
                    np.log(1e3), np.log(fmax), epsabs=0, epsrel=1e-12,
                    limit=500)[0]
    i1, i2 = names.index('o1'), names.index('o2')
    assert np.max(np.abs(K['mixed'] - K['split'])) <= \
        1e-12 * np.max(np.abs(K['split'])), 'the split is not exact'
    for (i, j), ef, tol in (((i1, i1), 0.8, 1e-9), ((i2, i2), 2.0, 1e-9)):
        assert abs(K['mixed'][i, j] / ref(ef) - 1.0) < tol, (ef, K['mixed'][i, j])
    for (i, j), ef in (((i1, i1), 0.8), ((i2, i2), 2.0), ((i1, i2), 1.4)):
        assert abs(abs(K['corr'][i, j]) / ref(ef) - 1.0) < 1e-6, \
            (ef, K['corr'][i, j] / ref(ef))


def test_a_coloured_covariance_integrates_any_flicker_exponent():
    """The band integral for a flicker exponent other than 1 (2026-09-25).
    `_coloured_covariance` took the trapezoid in ln nu, exact for EF = 1
    only: ``((1 - EF) ln r)^2 / 12`` off on a flat response, +2.78e-4 at
    EF = 2 and 40 per decade (predicted 2.76e-4).  The product rule (the
    power law analytic per interval) removes that, and leaves the response's
    bend at ``h^2`` -- a term that telescopes at EF = 1 but not here
    (+1.71e-5 / -2.70e-6 at EF = 0.8 / 2.0, 4.0x per halving of the grid
    ratio); Richardson on the same points removes that too.  Measured, RC
    driven with `flicker_noise(k, EF)` against the closed-form integrand
    by adaptive quadrature: -6.1e-10 (EF = 0.8) and +1.2e-9 (EF = 2.0) at
    40 per decade, fourth order down to the method's floor."""
    from scipy.integrate import quad
    import warnings
    circuit.default_toolkit = circuit.numeric
    T, Rv, Cv, k = 1e-6, 1e3, 1e-9, 1e-20
    fc = 1.0 / (2.0 * np.pi * Rv * Cv)
    for ef in (0.8, 2.0):
        c = SubCircuit()
        c.add_node('in')
        c.add_node('out')
        c['V'] = VSin('in', gnd, va=0.1, vo=0.0, freq=1.0 / T)
        c['R'] = R('in', 'out', r=Rv, noisy=False)
        c['C'] = C('out', gnd, c=Cv)
        c['n'] = _pow_flicker(ef)('out', gnd, k=k)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        o = [str(x) for x in c.nodes if str(x) != 'gnd!'].index('out')
        fmin, fmax = 1e3, 0.5 * len(pss.factored_period().steps) / T
        K = pac.covariance(pss, colour_fmin=fmin)
        ref = quad(lambda u: k * np.exp((1.0 - ef) * u) * Rv ** 2
                   / (1.0 + np.exp(2.0 * u) / fc ** 2),
                   np.log(fmin), np.log(fmax), epsabs=0, epsrel=1e-12,
                   limit=500)[0]
        assert abs(K[o, o] / ref - 1.0) < 1e-8, (ef, K[o, o] / ref - 1.0)


class _ModLorentz(IS):
    """A Lorentzian current whose level follows the terminal voltage --
    coloured, not a power law, AND modulated: refused by the band
    integral."""

    def CY(self, x, w, epar=None):
        v = float(x[0] - x[1])
        p = self.iparv.noisePSD * (1.0 + v * v) / (1.0 + (float(w) * 1e-7) ** 2)
        return self.toolkit.array([[p, -p], [-p, p]])


def test_a_coloured_covariance_integrates_a_stationary_lorentzian_source():
    """A colour that is NOT a power law (2026-09-25): `IS(noiseTau)`, a
    Lorentzian ``P / (1 + (w tau)^2)``.  `NoiseComponents.model` calls it
    per-band, and the band integral refused it; a STATIONARY one enters
    through its own `CY(nu)` and unit sources on its support (no square
    root: `_coloured_covariance`).

    Two references, measured (radau):
      * a driven RC, against the EXACT band-limited closed form
        ``int P R^2 / ((1 + (2 pi f tau)^2)(1 + (2 pi f RC)^2)) df``:
        -1.2e-7 at 100 / 200 / 400 points; and against the SAME noise
        realised as a white source through an explicit RC filter into a
        transconductor -- the white Lyapunov route, 3e-15 .. 1.6e-11 of the
        full-band closed form -- -3.3e-6, exactly the band left below
        `fmin = 1e-6 f0` (predicted 3.2e-6).  Gear: both routes second
        order (-9.2e-5 / -1.6e-4 at 200 points).
      * the comparator-jitter sampler with the Lorentzian on its threshold
        node: the crossing's coloured sigma^2 against ``Var_band(v_n) /
        s_1^2`` -- the crossings' path (`dtheta`) of the same integral:
        +1.4e-5 at 100 points (tau = 0.3 T), +6.8e-7 at 200.  With the
        corner at 8 f0 (tau = 0.02 T, the node's own RC one step) it read
        +2.8e-4 / +1.7e-5: radau's transfer at that grid, not the integral.
    A Lorentzian MODULATED by the orbit was refused here until 2026-09-26;
    it is now integrated (the separable path)."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 1e-6

    def band(P, R_, RC, tau, f1, f2):
        a, b = 2.0 * np.pi * tau, 2.0 * np.pi * RC
        F = lambda f: (a * np.arctan(a * f) - b * np.arctan(b * f)) / (a * a - b * b)  # noqa: E731
        return P * R_ ** 2 * (F(f2) - F(f1))

    Rv, Cv = 1e3, 0.5e-9
    Rf, Cf, g, Pw = 1e3, 0.3e-9, 1e-3, 1e-20
    tau, P = Rf * Cf, g * g * Pw * Rf * Rf
    out = {}
    for kind in ('coloured', 'filtered'):
        c = SubCircuit()
        c.add_node('in')
        c.add_node('out')
        c['V'] = VSin('in', gnd, va=0.1, vo=0.0, freq=1.0 / T)
        c['R'] = R('in', 'out', r=Rv, noisy=False)
        c['C'] = C('out', gnd, c=Cv)
        if kind == 'coloured':
            c['n'] = IS('out', gnd, i=0.0, noisePSD=P, noiseTau=tau)
        else:
            c.add_node('f')
            c['nw'] = IS('f', gnd, i=0.0, noisePSD=Pw)
            c['rf'] = R('f', gnd, r=Rf, noisy=False)
            c['cf'] = C('f', gnd, c=Cf)
            c['gm'] = BSource('f', gnd, gnd, 'out', i_func=lambda u: g * u)
        pss = PSS(c, method='radau', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 200, maxiterations=40)
        pac = PAC(c, toolkit=circuit.numeric)
        o = [str(x) for x in c.nodes if str(x) != 'gnd!'].index('out')
        fN = 0.5 * len(pss.factored_period().steps) / T
        band_ = dict(colour_fmin=1e-6 / T) if kind == 'coloured' else {}
        out[kind] = (pac.covariance(pss, **band_)[o, o], fN)
    kc, fN = out['coloured']
    kf, _ = out['filtered']
    assert abs(kc / band(P, Rv, Rv * Cv, tau, 1e-6 / T, fN) - 1.0) < 1e-6
    full = band(P, Rv, Rv * Cv, tau, 0.0, np.inf)
    assert abs(kf / full - 1.0) < 1e-9
    assert abs(kc / kf - band(P, Rv, Rv * Cv, tau, 1e-6 / T, fN) / full) < 1e-6

    ## the crossings' path: the Lorentzian on the jitter sampler's threshold
    Rn, Cn, tn, Pn = 1e4, 1e-12, 3e-7, 5e-23
    s1 = 5.0 / (0.9 * T)
    sig = {}
    for flick in (False, True):
        cir = _jitter_sampler(T)
        if flick:
            cir['Fn'] = IS('n', gnd, i=0.0, noisePSD=Pn, noiseTau=tn)
        pss = PSS(cir, method='radau', reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 100, maxiterations=100,
                      state_events=True)
            band_ = dict(colour_fmin=1e-6 / T) if flick else {}
            sig[flick] = float(PAC(cir, toolkit=circuit.numeric).event_jitter(
                pss, **band_)['sigma_t'][0]) ** 2
        fN = 0.5 * len(pss.factored_period().steps) / T
    var_n = band(Pn, Rn, Rn * Cn, tn, 1e-6 / T, fN)
    assert abs((sig[True] - sig[False]) / (var_n / s1 ** 2) - 1.0) < 5e-5, \
        (sig, var_n / s1 ** 2)

    ## modulated AND not a power law: integrated since 2026-09-26
    c = SubCircuit()
    c.add_node('in')
    c.add_node('out')
    c['V'] = VSin('in', gnd, va=0.5, vo=1.0, freq=1.0 / T)
    c['R'] = R('in', 'out', r=Rv, noisy=False)
    c['C'] = C('out', gnd, c=Cv)
    c['n'] = _ModLorentz('out', gnd, i=0.0, noisePSD=P)
    pss = PSS(c, method='radau', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 100, maxiterations=40)
        ## (refused until 2026-09-26; a level under a fixed shape is now the
        ## separable path -- see
        ## `test_a_coloured_covariance_takes_a_modulated_non_power_law_source`)
        with pytest.warns(RuntimeWarning, match='SIGN-BLIND'):
            Km = PAC(c, toolkit=circuit.numeric).covariance(pss, colour_fmin=1e-6 / T)
        assert np.all(np.isfinite(Km)) and np.max(np.diag(Km)) > 0.0


def test_the_bordered_consumers_run_on_a_staged_gear_solve_too():
    """Item 2 of the 2026-09-22 list, finished (2026-09-23): on a staged
    GEAR solve the adjoint row and the covariance closure were unbordered
    -- the columns did not exist until item 3 stored them in gear's pair
    form, and both consumers assumed the one-step width.  Now: the adjoint
    row does the same elimination on the pair map with the event rows'
    term as a second injected reverse pass, and `_event_closure` builds
    its noise-sensitivity rows at the MAP's width (the pair's `(x_j,
    x_{j-1})`).

    Measured.  (a) The adjoint row on the PWM loop at 100 points is
    dual-consistent with the bordered forward solve to 1.7e-14 at l = 0
    and 7e-14 at l = 1, where the UNBORDERED row is 15 % / 10 % off --
    the same order the bordered PAC removed.  (b) The covariance closure
    on `_jitter_sampler`: gear's held variance tracks its OWN threshold
    variance to 3 % / 0.4 % / 0.2 % at 100 / 400 / 800 points.  ⚠ The
    ABSOLUTE numbers there are 0.60 / 0.85 / 0.92 of the analytic, and
    that is gear's known covariance floor, not the bordering: `Var(n)` is
    the threshold node's own kT/C with no event in it and reads the same
    0.60 / 0.85 / 0.92 whether the solve is staged or not (its RC is
    10 ns against a 10 ns step at N = 100; radau reads 1.0004 throughout).
    So the pin is the RATIO the bordering is responsible for.  ⚠ AND THE
    BORDERING IS NOT WHAT SAVES THE HELD VARIANCE HERE: with the grid
    landed but the columns withheld gear reads 0.607 against the bordered
    0.6016 -- 1 % -- while its SIDEBAND response is 10 % off unbordered
    (above).  What the hold needs is the landed grid: an UNSTAGED gear
    solve of the same circuit reads 537x the analytic held variance, and
    that is the contrast pinned below."""
    import warnings as _w
    from pycircuit.circuit import remove_row_col
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    f0 = 1.0 / T
    ## (a) the adjoint row on a staged gear solve
    cir = _pwm_loop(T)
    p = PSS(cir, method='gear', reltol=1e-9)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T, timestep=T / 100, maxiterations=100, state_events=True)
    assert p.converged and p._event_columns is not None
    names = [str(n_) for n_ in cir.nodes]
    io_full = names.index('out')
    io = io_full if io_full < p.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), p.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    fin = 0.3 * f0
    pac = PAC(cir, toolkit=circuit.numeric)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        res = pac.solve(p, freqs=[fin])
        H = np.asarray(pac.adjoint_sideband_row(p, fin, io, sidebands=[0, 1]))
        ev = p._event_columns
        p._event_columns = None
        try:
            H0 = np.asarray(pac.adjoint_sideband_row(p, fin, io, sidebands=[0, 1]))
        finally:
            p._event_columns = ev
    fout = np.asarray(res.sweep_values, dtype=float)
    X = np.asarray(res.x)
    for li, l in enumerate((0, 1)):
        k = int(np.argmin(np.abs(fout - abs(fin + l * f0))))
        x = complex(X[io_full, k])
        assert abs(x - complex(H[li] @ u_ac)) < 1e-10 * abs(x), (l, x, H[li] @ u_ac)
        assert abs(x - complex(H0[li] @ u_ac)) > 1e-2 * abs(x), (l, x, H0[li] @ u_ac)
    ## (b) the covariance closure in pair form
    kT = 1.380649e-23 * 300.0
    exp_n = kT / 1e-12
    exp_h = 4.0 * exp_n
    got = {}
    for N in (100, 400):
        c3 = _jitter_sampler(T)
        ps = PSS(c3, method='gear', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            ps.solve(period=T, timestep=T / N, maxiterations=100, state_events=True)
        red3 = [str(n_) for i, n_ in enumerate(c3.nodes) if i != ps.irefnode]
        ih, inn = red3.index('hold'), red3.index('n')
        ts3 = np.asarray(ps.waveform[0], dtype=float)
        j7 = int(np.searchsorted(ts3, 0.7 * T))
        pac3 = PAC(c3, toolkit=circuit.numeric)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            K0, seq = pac3.covariance(ps, samples=True)
            ## the contrast that matters to a user: NOT staging at all
            pu = PSS(_jitter_sampler(T), method='gear', reltol=1e-9)
            pu.solve(period=T, timestep=T / N, maxiterations=100, state_events=False)
            _K0u, sequ = PAC(pu.cir, toolkit=circuit.numeric).covariance(pu, samples=True)
        ## m x m whatever the method (2026-09-29); gear's pair on request
        assert np.shape(K0) == (c3.n - 1, c3.n - 1)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            assert np.shape(pac3.covariance(ps, pair=True)) == (2 * (c3.n - 1),) * 2
        var_n = float(np.mean([K[inn, inn] for K in seq])) / exp_n
        held = float(seq[j7][ih, ih]) / exp_h
        got[N] = held / var_n
        ## the bordering's job is the ratio; the level is gear's own floor
        assert 0.5 < var_n < 1.0, (N, var_n)
        assert abs(got[N] - 1.0) < 5e-2, (N, got[N], var_n, held)
        ## an UNSTAGED gear solve of the same circuit: nonsense at the hold
        tsu = np.asarray(pu.waveform[0], dtype=float)
        assert float(sequ[int(np.searchsorted(tsu, 0.7 * T))][ih, ih]) / exp_h > 100.0
    assert abs(got[400] - 1.0) < abs(got[100] - 1.0), got


@pytest.mark.parametrize('method', ['radau', 'gear', 'trap', 'glm2'])
def test_the_state_event_stage_runs_matrix_free(method):
    """The state-event stage under `matrix_free=True` (2026-09-25; until then
    it warned and was skipped, the crossings left inside their steps).  Its
    bordered system is a mat-vec: one forward replay per Krylov direction
    gives the map and the state at every event node; the walk is factored
    and carries only the event (and period) columns; `_monodromy` stays
    None; the event rows' derivatives `G` come from reverse replays and the
    dense map to every node (`P_nodes`) is built on first read.

    Three things had to change in the matrix-free Newton to get there, each
    measured:
    * a LINE SEARCH, as the dense path's `fsolve(line_search=True)`: the
      UNSTAGED matrix-free solve of the PWM loop failed under every method,
      radau included, and a first staged step moved two crossings by a
      whole period;
    * COLUMN SCALING of the crossing and period columns (``dx/dT`` ~1e6):
      the comparator oscillator's bordered operator had condition 3.7e6,
      16 scaled;
    * an inner solve judged by the RESIDUAL it reached: glm2's staged
      system (scaled condition 1.6e6) stalled at 4.4e-11 against a 1e-11
      target -- a failure flag on an excellent step.

    Measured against the dense staged solves: the crossings to 1.6e-15 ..
    9.7e-13, the periods to 3.6e-14, the waveforms to 1.4e-10 .. 3.7e-10
    (PWM) and 1e-14 .. 9e-13 (oscillator).  Under radau the consumers too:
    `G` 7.7e-15, PAC 3.9e-13, the adjoint row 1.8e-12, the covariance (its
    closure reads `P_nodes`) 1.2e-14; the PPV samples 2.8e-11 and the
    multipliers 3e-15 on the oscillator.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    Tp = 1e-5
    f0 = 1.0 / Tp

    def pwm():
        cir = _pwm_loop(Tp)
        del cir['Vin']
        cir.add_node('vin0')
        cir['Vin'] = VS('vin0', gnd, v=5.0, vac=0.0)
        cir['Vp'] = VSin('vin', 'vin0', vo=0.0, va=0.0, freq=0.3 * f0, phase=0.0, vac=1.0)
        cir['Vramp'].iparv.vac = 0.0
        return cir

    def rel(a, b):
        a, b = np.asarray(a), np.asarray(b)
        return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))

    seed, Tl = _relaxation_oscillator_seed(_comparator_relaxation_oscillator())
    cases = ((pwm, dict(period=Tp, timestep=Tp / 60, x0=np.zeros(pwm().n - 1))),
             (_comparator_relaxation_oscillator, dict(period=Tl, timestep=Tl / 200, x0=seed)))
    for mk, kw in cases:
        ## ⚠ trap's PWM loop SETTLED FIRST (`tstab`) since a landed ramp edge
        ## drops the order (2026-09-28): from zeros its staged Newton stalls
        ## at |F| 3.8 (1017 evaluations; the map is smooth -- the same steps
        ## drop in every walk -- and its event columns FD-exact), where
        ## from a settled seed it converges in 14 to the same crossings
        ## (`_staged_fallback` recovers the zero start since, at ~9x the time)
        if mk is pwm and method == 'trap':
            kw = dict(kw, tstab=20 * Tp)
        got = {}
        for mf in (False, True):
            cir = mk()
            p = PSS(cir, method=method, reltol=1e-10 if mk is pwm else 1e-9)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                p.solve(maxiterations=100, matrix_free=mf, **kw)
            assert p.converged and p._event_columns is not None, (method, mk.__name__, mf)
            got[mf] = (p, cir)
        (pd, cd), (pm, cm) = got[False], got[True]
        assert pm._monodromy is None and pd._monodromy is not None
        assert rel(pm._state_event_fracs, pd._state_event_fracs) < 1e-10, method
        assert abs(pm.period / pd.period - 1.0) < 1e-12, method
        assert rel(pm.waveform[1], pd.waveform[1]) < 1e-8, method
        assert rel(pm._event_columns['G'], pd._event_columns['G']) < 1e-9, method
        if method != 'radau':
            continue
        if mk is pwm:
            out = {}
            for mf, (q, c) in got.items():
                pac = PAC(c, toolkit=circuit.numeric)
                with _w.catch_warnings():
                    _w.simplefilter('ignore')
                    out[mf] = (np.asarray(pac.solve(q, [0.3 * f0]).x),
                               pac.adjoint_sideband_row(q, 0.3 * f0, 1, sidebands=[0, 1]),
                               pac.covariance(q))
            for a, b in zip(out[True], out[False]):
                assert rel(a, b) < 1e-9, method
        else:
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                sd = np.asarray(pd.ppv()[1]['samples'])
                sm_ = np.asarray(pm.ppv()[1]['samples'])
            assert rel(sm_, sd) < 1e-8, method
    ## the unstaged matrix-free solve of the PWM loop (it failed under every
    ## method before the line search)
    cir = pwm()
    p = PSS(cir, method=method, reltol=1e-9)
    q = PSS(pwm(), method=method, reltol=1e-9)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=Tp, timestep=Tp / 100, x0=np.zeros(cir.n - 1), maxiterations=100,
                matrix_free=True, state_events=False)
        q.solve(period=Tp, timestep=Tp / 100, x0=np.zeros(cir.n - 1), maxiterations=100,
                state_events=False)
    assert p.converged, method
    assert rel(p.waveform[1], q.waveform[1]) < 1e-10, method


def test_a_glm_run_reads_its_small_signal_off_its_own_map_and_its_covariance_off_a_twin():
    """What a Nordsieck GLM run's state-space consumers read (2026-09-25,
    Andreas: "Native for all but covariance"):

    * DRIVEN, `PAC.solve`, the adjoint rows, pnoise, `sampled_noise`: the
      GLM's own map on the state (`PSS._state_map`) -- exact at its order
      (`test_a_glm_reads_pac_off_its_own_map_at_its_order`), so no longer
      equal to a radau solve's.
    * The covariance surfaces: a radau twin by default (`_lyapunov_host`),
      equal to a radau solve's on the same grid; the GLM's own with
      `monodromy='native'` (`PAC._lyapunov_pieces_glm`: one shared sample
      per step, first order -- the effective stage weights ``l^T B`` are
      not all positive).  The sampler's held variance: glm2 1.0153 /
      1.0055 / 1.0022 kT/C at 200 / 400 / 800 points (radau twin
      1.000000).
    * An OSCILLATOR's small-signal surfaces read the GLM's own map too
      (Andreas: "If GLM is accurate use GLM"), with the deflated answer
      UNREFINED (`PAC._deflated_solve`), MEASURED: the GLM map's unit
      multiplier sits ``eta = O(h^p)`` off 1 (its startup breaks the phase
      symmetry; van der Pol in LC form with a ``0.3 u^2`` asymmetry, glm3
      at 60 points: 2.35e-5), and refined on that operator the answer near
      a harmonic missed the physical one by ``eta / (2 pi r)`` -- 3.8e-3 /
      0.35 / 1.0 at r = 1e-3 / 1e-5 / 1e-7, the pole lost.  Unrefined, the
      carrier's response against radau at 480 points: 3.3e-5 / 3.3e-5 /
      8.4e-5, forward and adjoint alike (dual to 1e-6, O(eta)).  ⚠ A weak
      sideband is swamped: the l = -1 one here is 2e6 below the carrier's,
      and the GLM's error there, 2e-6 of the carrier's response at its
      order, is several times the coefficient (radau: 1e-4).
    """
    import warnings as _w
    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric

    def build():
        c = _q20_rlc()
        c['n'] = IS('c', gnd, i=0.0, noisePSD=1e-20)
        return c

    def run(method, mono=None):
        c = build()
        p = PSS(c, method=method, reltol=1e-10)
        if mono is not None:
            p.monodromy = mono
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=1e-3, timestep=1e-3 / 100)
            pac = PAC(c, toolkit=circuit.numeric)
            res = pac.solve(p, [300.0])
            K0 = np.asarray(pac.covariance(p), dtype=float)
        k = int(np.argmin(np.abs(np.asarray(res.sweep_values, dtype=float) - 300.0)))
        x = complex(np.asarray(res.x)[[str(n_) for n_ in c.nodes].index('c'), k])
        return x, K0, p

    x_r, K_r, _pr = run('radau')
    for method in ('glm2', 'glm3'):
        x, K, p = run(method)
        ## PAC is the GLM's own: not radau's, and at its order (300 Hz:
        ## glm2 2.6e-6, glm3 5.3e-9 against AC at 100 points)
        assert 1e-13 * abs(x_r) < abs(x - x_r) < 1e-5 * abs(x_r), (method, x, x_r)
        ## the covariance is the radau twin's
        assert np.max(np.abs(K - K_r)) < 1e-12 * np.max(np.abs(K_r)), method
        assert p.factored_period().kind == 'glm' and p.monodromy_twin() is p
        ## native: the GLM's own Lyapunov pieces, first order (glm2 0.164
        ## off radau's at 100 points on this Q = 20 resonator)
        x_n, K_n, _pn = run(method, 'native')
        assert x_n == x
        assert 0.02 < np.max(np.abs(K_n - K_r)) / np.max(np.abs(K_r)) < 0.4, method

    ## the sampler: the native covariance's held value converges to kT/C
    ktc = _KB * _TEMP / 100e-12
    held = []
    for npts in (200, 400):
        cir, pss, io, pac, T = _sampler_fixture_method(
            lambda c: c.__setitem__('S0', _sw()), 'glm2', npts)
        N = len(pss.factored_period().steps)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            K0t, seq_t = pac.covariance(pss, samples=True)
            pss.monodromy = 'native'
            K0n, seq_n = pac.covariance(pss, samples=True)
        jh = int(0.375 * N)
        assert abs(seq_t[jh][io, io] / ktc - 1.0) < 1e-4
        held.append(seq_n[jh][io, io] / ktc - 1.0)
    assert 0.0 < held[1] < held[0] < 0.03 and held[0] / held[1] > 2.0, held

    ## an oscillator's small-signal surfaces: the GLM's own map, the pole
    ## carried analytically (an asymmetric van der Pol in LC form)
    def vdp():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: u - u ** 3 / 3.0 + 0.3 * u ** 2)
        c['ac'] = IS('v', gnd, i=0.0, iac=1.0)
        return c

    def carrier_response(method, N):
        cir = vdp()
        q = PSS(cir, method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=6.66, timestep=6.66 / N, x0=np.array([2.0, 0.0]),
                    maxiterations=100)
        assert q.converged
        f0 = 1.0 / float(q.period)
        iv = [str(n_) for n_ in cir.nodes].index('v')
        io = iv if iv < q.irefnode else iv - 1
        (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), q.irefnode, circuit.numeric)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        out = []
        for r in (1e-3, 1e-5, 1e-7):
            f = (1.0 + r) * f0
            pac = PAC(cir, toolkit=circuit.numeric)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                res = pac.solve(q, [f], sweeptype='absolute')
                h = complex(np.asarray(pac.adjoint_sideband_row(q, f, io, sidebands=[0]))[0] @ u_ac)
            sv = np.asarray(res.sweep_values, dtype=float)
            X = np.asarray(res.x)[iv]
            ks = [k for k in range(len(sv)) if abs(sv[k] - f) < 1e-9 * f]
            out.append((complex(X[max(ks, key=lambda k: abs(X[k]))]), h))
        return q, out

    q, got = carrier_response('glm3', 60)
    assert q._state_map().is_glm and q.monodromy_twin() is q
    _qr, ref = carrier_response('radau', 480)
    for (x, h), (xr, hr), r in zip(got, ref, (1e-3, 1e-5, 1e-7)):
        assert abs(x / xr - 1.0) < 2e-4, (r, x, xr)
        assert abs(h / hr - 1.0) < 2e-4, (r, h, hr)
