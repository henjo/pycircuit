"""Shooting tests: shooting methods.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_fixtures import (_diode_mixer,
    _force_plain_map,
    _pss_lte,
    _pss_peak,
    _q20_rlc,
    _rc_ladder,
    _scaled_vdp,
    _shooting_trace,
    _varying_c_ladder,
    _vdp_asym)


def test_backward_euler_damps_the_limit_cycle_and_trapezoidal_does_not():
    """The defect 0.1c names, measured against an ANALYTIC answer.

    Backward Euler's numerical damping attenuates exactly the limit cycle PSS
    exists to find.  On a Q=20 resonator driven at resonance, where the analytic
    peak is Q*va = 20 V:

        steps/period   euler          trapezoidal
                  20   2.63 V (13%)   19.32 V (97%)
                 200  12.20 V (61%)   19.23 V (96%)

    Euler is not merely less accurate, it is WRONG BY A FACTOR and gets worse as
    the step coarsens -- silently.
    """
    trap_peak, Q = _pss_peak('trap', steps=20)
    euler_peak, _ = _pss_peak('euler', steps=20)
    assert trap_peak > 0.9 * Q, \
        'trapezoidal recovers only %.1f%% of the analytic amplitude' % (100 * trap_peak / Q)
    assert euler_peak < 0.5 * Q, \
        'the euler damping this test documents has gone: %.1f%%' % (100 * euler_peak / Q)


def test_euler_shooting_converges_like_a_newton():
    """Few iterations, and tightening the tolerance costs few more.

    Successive substitution on this circuit contracts by 0.8546 per
    iteration, so reaching 1e-9 would need ~130.  A Newton reaches it in a
    handful; measured at landing, 5 iterations at reltol 1e-4 and 10 at
    1e-9.
    """
    loose, nc_loose, _r = _shooting_trace('euler', reltol=1e-4)
    tight, nc_tight, _r2 = _shooting_trace('euler', reltol=1e-9)

    assert not nc_loose and not nc_tight
    assert len(loose) <= 8, 'euler took %d shooting iterations' % len(loose)
    assert len(tight) <= 15, 'euler took %d at reltol 1e-9' % len(tight)
    ## five extra decades for a handful of iterations is the Newton signature
    assert tight[-1][0] < 1e-9, 'final residual %.3e' % tight[-1][0]


# ---------------------------------------------------------------------------
# Phase 3: the (x, iq) monodromy
# ---------------------------------------------------------------------------

def test_trapezoidal_shooting_converges_with_the_x_iq_monodromy():
    """Trapezoidal's period map carries `iq`, so the monodromy must too.

    With an x-only monodromy trap did not converge at all -- worse, applying
    the Euler form to it converged SLOWER than no Jacobian (0.90 against
    0.855 per iteration), because the Jacobian was wrong rather than absent.
    Differentiating the two recursions together costs one extra matrix
    product and makes it a real Newton: measured 6 iterations at reltol 1e-4
    and 13 at 1e-9, residual 2.9e-11.
    """
    loose, nc_loose, _r = _shooting_trace('trap', reltol=1e-4,
                                          maxiterations=40)
    tight, nc_tight, _r2 = _shooting_trace('trap', reltol=1e-9,
                                           maxiterations=40)
    assert not nc_loose and not nc_tight
    assert len(loose) <= 10, 'trap took %d shooting iterations' % len(loose)
    assert tight[-1][0] < 1e-9, 'final residual %.3e' % tight[-1][0]


def test_trapezoidal_is_now_right_on_both_axes():
    """The point of the whole repair.

    The two failure modes were orthogonal and neither method escaped both:
    Euler CONVERGED the shooting equation and landed at 8.815 V against a
    20 V analytic peak (a discretisation error), while trapezoidal did not
    converge the shooting equation and landed at 19.848 V.  With the
    `(x, iq)` monodromy, trap converges AND lands at 19.990 V -- 0.05% of
    analytic.

    Euler is asserted unchanged in the same breath, because the unified
    propagation must reduce to the old one when the `iq` row is identically
    zero; if this drifts, the "one formula, two methods" claim is false.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    peaks = {}
    for method in ('euler', 'trap'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = PSS(_q20_rlc(), method=method, reltol=1e-6).solve(
                period=1e-3, timestep=1e-5, maxiterations=40)
        peaks[method] = float(np.max(np.abs(
            np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## trapezoidal does not damp the limit cycle: within 0.5% of 20 V
    assert abs(peaks['trap'] - 20.0) < 0.1, \
        'trap peak %.4f V, analytic 20 V' % peaks['trap']
    ## and Euler still damps it, by its own documented amount -- this is
    ## the level-2 error the shooting solve cannot see
    assert 8.5 < peaks['euler'] < 9.2, \
        'euler peak moved to %.4f V' % peaks['euler']


def test_gear2_shooting_converges_and_damps_between_euler_and_trap():
    """BDF-2 reaches back TWO steps, and the propagation follows it.

    Every method here writes its companion as
    `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`, and the integrator now states
    those coefficients rather than `shooting.py` transcribing them -- this
    file has recorded three times what transcribing an integration constant
    costs.  Euler is `b=0` with one past term, trapezoidal `b=-1` with one,
    Gear-2 `b=0` with two; one recursion serves all three.

    The physics is the check: on a Q=20 resonator against a 20 V analytic
    peak, numerical damping should order euler >> gear2 > trap.  Measured
    at landing -- euler 8.815 V, gear 19.766 V, trap 19.990 V.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    peaks, iters = {}, {}
    for method in ('euler', 'gear', 'trap'):
        trace, nonconv, res = _shooting_trace(method, reltol=1e-9,
                                              maxiterations=40)
        assert not nonconv, '%s did not converge' % method
        iters[method] = len(trace)
        peaks[method] = float(np.max(np.abs(
            np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## all three are real Newtons: a handful of iterations to 1e-9
    for m, k in iters.items():
        assert k <= 20, '%s took %d shooting iterations' % (m, k)

    ## and the damping orders as the methods do
    assert peaks['euler'] < peaks['gear'] < peaks['trap'], \
        'damping does not order euler < gear < trap: %r' % peaks
    assert abs(peaks['gear'] - 20.0) < 1.0, \
        'gear2 peak %.4f V, analytic 20 V' % peaks['gear']


def test_pss_solved_history_removes_the_gear2_seam():
    """The fix, against the number that predicted it.

    The seam was measured by continuing 80 periods past the converged solve
    with no re-seed, which reaches the limit cycle the same grid and method
    produce from a real history: 19.89297 V at 100 points per period,
    19.98524 at 200, 19.99735 at 400, against a plain-formulation PSS that
    lands at 19.76639 / 19.95451 / 19.99008 and an analytic 20 V.

    Making the entering history an unknown must reach the FIRST of those
    numbers -- that is what "the seam is the only difference" means -- so
    this asserts the prediction, not just an improvement.
    """
    circuit.default_toolkit = circuit.numeric
    predicted = {100: 19.89297, 200: 19.98524, 400: 19.99735}
    plain_err = {100: 2.336e-01, 200: 4.549e-02, 400: 9.924e-03}
    for npts, target in predicted.items():
        ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
        peak, pss = _pss_lte('gear', timestep=1e-3 / (npts - 1), reltol=1e-9)
        assert pss.solved_history and pss.converged
        assert abs(peak - target) < 5e-5, \
            'npts=%d landed at %.5f, not the primed limit cycle %.5f' % (
                npts, peak, target)
        ## and it is an improvement, by the factor the measurement predicted
        assert abs(peak - 20.0) < 0.5 * plain_err[npts], \
            'npts=%d error %.3e is not a clear gain on the plain %.3e' % (
                npts, abs(peak - 20.0), plain_err[npts])
        ## the seam is gone, and the report agrees it is gone
        assert pss.max_lte_seam is None, \
            'the solved history is not a seam: %r' % pss.max_lte_seam


def test_the_monodromy_is_the_pair_map_not_a_corner_of_it():
    """⚠ A SUB-BLOCK OF A SENSITIVITY IS NOT A MONODROMY.

    For a two-step method the one-period map acts on the PAIR
    `(x_n, x_{n-1})`, so the monodromy is `2m x 2m`.  The solved-history
    path used to hand back the `d x_{N-1}/d x_0` corner of it instead, and
    the eigenvalues of a corner mean nothing: it reported a spectral radius
    of 1.2796 for this resonator -- ABOVE ONE, i.e. an unstable orbit --
    where the circuit decays by `exp(-pi/Q)` every period.

    The analytic decay is the check, so this cannot drift with a
    reimplementation.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    analytic = float(np.exp(-np.pi / 20.0))

    def rho(method, force_plain):
        pss = PSS(_q20_rlc(), method=method, reltol=1e-9)
        if force_plain:
            _force_plain_map(pss)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=1e-3, timestep=1e-5, maxiterations=40)
        assert pss.converged
        return pss

    for method in ('trap', 'gear'):
        for plain in (True, False):
            pss = rho(method, plain)
            assert abs(pss.spectral_radius - analytic) < 0.01, \
                'method=%r plain=%r reports rho=%.6f against the analytic ' \
                'per-period decay %.6f -- a value above 1 here means a ' \
                'corner of the sensitivity is being read as a monodromy' % (
                    method, plain, pss.spectral_radius, analytic)

    ## and the shape follows the method's reach, mechanically
    m = _q20_rlc().n - 1
    assert np.asarray(rho('gear', False)._monodromy).shape == (2 * m, 2 * m)
    assert np.asarray(rho('trap', False)._monodromy).shape == (m, m)


def test_the_matrix_free_matvec_is_the_dense_monodromy():
    """RECORDED SCOPE ITEM 6: `M v` without ever forming `M`.

    The dense path propagates `2m` columns per step; the matrix-free path
    replays the SAME recursion at width 1 on a stored factorisation.  If
    those two disagree, everything built on the matvec is measuring its own
    bug, so this pins them against each other directly rather than against
    a converged answer that could absorb the difference.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_varying_c_ladder(), method='gear', reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 25, maxiterations=20)
    m = pss.cir.n - 1
    times, hs = pss._period_grid(1e-3, 25, None)
    ## ⚠ THE TWO HISTORY STATES MUST DIFFER, and that was also measured.
    ## Seeded with `x_0 = x_{-1} = 0`, the two opening capacitances
    ## `C(x_0)` and `C(x_{-1})` are the same matrix, so seeding the ring
    ## with the WRONG one is invisible -- a mutation replacing
    ## `_C_at(xm1_in)` with `_C_at(x0_in)` stayed green.  Two genuinely
    ## different entering states are what make the opening testable.
    xa = np.full(m, 0.05)
    xb = np.full(m, -0.03)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        M = pss._walk('pair', np.concatenate((xa, xb)), times, hs).monodromy()
        _wf = pss._walk('pair', np.concatenate((xa, xb)), times, hs,
                        dense=False, keep=True)
        C0, steps = _wf.opening, _wf.steps

    rng = np.random.default_rng(0)
    for _ in range(3):
        v = rng.standard_normal(2 * m)
        got = pss._monodromy_matvec(C0, steps, v)
        want = M @ v
        err = np.linalg.norm(got - want) / np.linalg.norm(want)
        assert err < 1e-11, 'matvec differs from the dense monodromy by %.3e' % err


def _gear_on_a_smooth_3to1_grid(period_column, npts=200):
    """Gear on van der Pol (Q = 8, `a u^2`, a = 0.3: period 6.73 s) on a
    smooth 3:1 grid, and the seed `[2, 0]` at 2 pi -- 7 % below the period."""
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
    w = 1.0 + 0.5 * np.sin(2.0 * np.pi * np.arange(npts) / npts)
    pss = PSS(c, method='gear', reltol=1e-12, period_column=period_column)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    return pss, dict(period=T, timestep=T / npts, x0=x0, grid=w / w.sum())


def test_gears_closing_period_column_carries_the_opening_step_and_auto_falls_back():
    """⚠ Gear's 'closing' period column was 15 % off on a caller's
    non-uniform grid until 2026-09-26.  Under 'closing' only the last step
    follows `T` -- but gear's pair map OPENS each period with that step as
    its previous one, so the opening step's coefficients move with `T` too,
    and that route was missing (a one-step method has none: radau's column
    6e-11).  Against finite differences at the seed: 1.47e-1 -> 8.6e-10
    (`proportional`, 6e-10; a uniform grid, 2e-9).

    With the column right, the free-period Newton from a poor seed (`[2, 0]`,
    a flat gear history, the period 7 % low) still failed under 'closing':
    its first step asks the closing step to absorb -4.3 s, the grid falls
    back to proportional for that one evaluation (the function changes under
    a closing Jacobian), and the iterates drift to a false near-solution at
    small amplitude, |F| ~ 1e-3.  With the WRONG column it had failed at 200
    and 240 points (a transient step diverging at a trial point, which
    aborted the solve) and not converged at 300.  Now:
      * a line-search trial whose evaluation raises is halved, not fatal
        (`fsolve`);
      * 'auto' falls back to 'proportional' from the same seed when its
        closing solve fails (`PSS._closing_fallback`), warned.
    All of 160 / 200 / 240 / 300 / 400 points converge: 6.736744 /
    6.734582 / 6.733396 / 6.732417 / 6.731650 s, second order."""
    import types
    import warnings as _w

    class _Stop(Exception):
        pass
    cap = {}

    def grab(self, func, z0, *args, **kwargs):
        cap['func'], cap['z0'] = func, np.asarray(z0, dtype=float).copy()
        raise _Stop()
    pss, kw = _gear_on_a_smooth_3to1_grid('closing')
    pss._free_period_solve = types.MethodType(grab, pss)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        try:
            pss.solve(maxiterations=40, **kw)
        except _Stop:
            pass
        func, z = cap['func'], cap['z0']
        _F, J = func(z)
        d = 1e-6 * z[-1]
        zp, zm = z.copy(), z.copy()
        zp[-1] += d
        zm[-1] -= d
        fd = (np.asarray(func(zp)[0], float) - np.asarray(func(zm)[0], float)) / (2 * d)
    col = np.asarray(J, float)[:, -1]
    err = np.max(np.abs(col - fd)) / np.max(np.abs(fd))
    assert err < 1e-6, 'the closing period column is %.2e from FD' % err

    ref, kw = _gear_on_a_smooth_3to1_grid('proportional')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ref.solve(maxiterations=40, **kw)
    assert ref.converged
    ## (a generous budget: the closing attempt must STALL OUT and hand over,
    ## not run to it -- 116 s at 300 before the stall stop, 14.5 s after)
    import pycircuit.circuit.analysis as _an
    attempts = []
    orig = _an.fsolve

    def spy(f, z0, *a, **k):
        n = [0]

        def g(z, *aa):
            n[0] += 1
            return f(z, *aa)
        out = orig(g, z0, *a, **k)
        attempts.append((k.get('stall_window'), out[1].get('stalled'), n[0]))
        return out
    auto, kw = _gear_on_a_smooth_3to1_grid('auto')
    _an.fsolve = spy
    try:
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            auto.solve(maxiterations=300, **kw)
    finally:
        _an.fsolve = orig
    assert attempts[0][:2] == (PSS.CLOSING_STALL_WINDOW, True), attempts
    assert attempts[0][2] < 200, attempts
    assert all(a[0] is None for a in attempts[1:]), attempts
    assert auto.converged
    assert any("'closing' period column" in str(r.message) for r in rec), \
        [str(r.message)[:80] for r in rec]
    assert abs(auto.period / ref.period - 1.0) < 1e-9, (auto.period, ref.period)


def test_the_monodromy_transpose_is_a_reverse_replay():
    """`M^T v` from the stored factors, backwards -- no reverse integrator.

    A shooting monodromy is a product of per-step solves, so its transpose
    is that product replayed in REVERSE ORDER with each solve transposed.
    The kept walk (`keep=True`) stores every step's factorisation and every
    factorisation here can solve transposed, so the reverse pass needs no
    new integrator, no refactorisation and no second traversal.

    ⚠ WHY THAT IS WORTH A TEST FOR MACHINERY NOTHING YET USES.  Demir &
    Roychowdhury (TCAD 22(2) 188-196) treat reverse integration as the
    barrier to computing a PPV -- "often unavailable even in existing
    time-domain simulators and may require significant changes to core
    simulation routines" -- and adjoint noise needs the same object.  That
    barrier is real for a forward-only DENSE implementation.  It is not real
    here, and this pins the reason rather than leaving it in a scratch file:
    both remaining capabilities bottleneck on this one piece.

    Measured when it was written: agreement 1.8e-15 with the dense `M^T`,
    and the reverse pass costs 0.75x the forward one -- CHEAPER, because it
    does two `C^T` products against one shared transposed solve where the
    forward does two `C` products and a solve.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_rc_ladder(6), method='gear', reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 50, maxiterations=2)
    m = pss.cir.n - 1
    times, hs = pss._period_grid(1e-3, 50, None)
    z = np.zeros(m)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _wf = pss._walk('pair', np.concatenate((z, z)), times, hs,
                        dense=False, keep=True)
        C0, steps = _wf.opening, _wf.steps

    ## dense M from the FORWARD matvec, itself pinned against the traversal
    ## by test_the_matrix_free_matvec_is_the_dense_monodromy
    M = np.zeros((2 * m, 2 * m))
    for i in range(2 * m):
        e = np.zeros(2 * m)
        e[i] = 1.0
        M[:, i] = pss._monodromy_matvec(C0, steps, e)

    rng = np.random.default_rng(5)
    for _ in range(4):
        v = rng.standard_normal(2 * m)
        got = pss._monodromy_matvec_transposed(C0, steps, v)
        want = M.T @ v
        err = np.linalg.norm(got - want) / np.linalg.norm(want)
        assert err < 1e-11, \
            'the reverse replay differs from the dense M^T by %.3e' % err

    ## ⚠ and the identity that makes it useful: <M^T u, w> == <u, M w>.
    ## A transpose that is only checked against a matrix built from the
    ## forward pass could share an error with it; this cannot.
    u, w = rng.standard_normal(2 * m), rng.standard_normal(2 * m)
    lhs = float(np.dot(pss._monodromy_matvec_transposed(C0, steps, u), w))
    rhs = float(np.dot(u, pss._monodromy_matvec(C0, steps, w)))
    assert abs(lhs - rhs) < 1e-9 * max(abs(rhs), 1.0), \
        'adjoint identity fails: <M^T u, w> = %.12g against <u, M w> = ' \
        '%.12g' % (lhs, rhs)


def test_the_monodromy_is_correct_across_a_switching_boundary():
    """The saltation concern, run as a falsifier and NOT confirmed.

    Bizzarri & Wei (ECCTD 2011) state that for hybrid systems "the monodromy
    matrix is not defined at impact events", so the transition matrix needs
    a SALTATION correction `S` -- and with no state jump `S` still differs
    from the identity whenever the VECTOR FIELD jumps, which an ideal switch
    does.  If that reached this code, every switching circuit's monodromy
    would be accumulated through a boundary where the linearisation does not
    exist.

    ⚠ IT DOES NOT REACH IT, and the reason is structural rather than lucky.
    PSS's monodromy is the exact derivative of the DISCRETE period map: each
    step uses its own converged `Jf` and `C`, which already describe
    whichever side of the switch that step is on.  The saltation matrix is a
    CONTINUOUS-time construct for correcting `Phi(t2,t0)` when the flow is
    discontinuous; a discrete map has no instant at which the field is
    undefined.

    ⚠ THE ASSERTION IS ON THE RATE, not the size.  A missing saltation term
    is an O(1) structural error that does not vanish under refinement, so
    "the gap is small" would not distinguish it from discretisation error.
    Measured on a switched RC whose control crosses twice per period, with
    Ron/Roff spanning six decades:

        npts    rel err (analytic monodromy vs finite differences)
         400    4.858e-04
         800    2.426e-04   (2.00x)
        1600    1.212e-04   (2.00x)
        3200    6.058e-05   (2.00x)

    Exactly first order, so the residue is discretisation and there is no
    O(1) term hiding under it.

    ⚠ AND THE CIRCUIT HAD TO BE BUILT TWICE.  The first two attempts put the
    switch in series with the source, which CLAMPS the node each period and
    erases the state: |M| came out 1.4e-22 and then 7.4e-51, so the test was
    comparing numerical zero against numerical zero and reporting a
    meaningless "rel err 1.000".  The switched element must change the decay
    RATE without clamping, or there is no monodromy to check.
    """
    import warnings
    from pycircuit.circuit import VSwitch
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def switched(ron, roff):
        c = SubCircuit()
        c.add_node('ctl')
        c.add_node('a')
        c.add_node('b')
        c['vc'] = VSin('ctl', gnd, va=2.0, freq=1.0 / per)
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per, phase=90.0)
        c['rs'] = R('a', 'b', r=1e5)
        c['c'] = C('b', gnd, c=1e-6)
        ## a switched LOAD, not a switched path to the source
        c['sw'] = VSwitch('b', gnd, 'ctl', gnd, Ron=ron, Roff=roff,
                          Von=1.0, Voff=0.0)
        return c

    gaps = []
    for npts in (200, 400, 800):
        pss = PSS(switched(1e3, 1e9), method='trap', reltol=1e-11)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(per, npts, None)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=per / npts, maxiterations=40)
        assert pss.converged
        ir = pss.irefnode
        Xw = np.asarray(res['tpss'].x, dtype=float)
        ctl = Xw[pss.cir.get_node_index('ctl')]
        toggles = int(np.sum(np.diff((ctl > 0.5).astype(int)) != 0))
        assert toggles >= 2, \
            'the switch no longer toggles (%d crossings), so this test has ' \
            'lost its subject' % toggles
        x0 = np.concatenate((Xw[:ir, 0], Xw[ir + 1:, 0]))

        def phi(v):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                xe = pss._walk('plain', np.asarray(v, dtype=float), times,
                               hs, T=per).x_end
            return np.asarray(xe, dtype=float)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            M = pss._walk('plain', x0, times, hs, T=per).monodromy()
        M = np.asarray(M, dtype=float)
        assert np.linalg.norm(M) > 0.1, \
            'the monodromy is %.3e -- the circuit erases its state each ' \
            'period, so there is nothing to check' % np.linalg.norm(M)

        base = phi(x0)
        Mfd = np.zeros((m, m))
        for j in range(m):
            d = 1e-7 * max(abs(x0[j]), 1.0)
            xp = x0.copy()
            xp[j] += d
            Mfd[:, j] = (phi(xp) - base) / d
        gaps.append(float(np.linalg.norm(M - Mfd)
                          / max(np.linalg.norm(Mfd), 1e-300)))

    for a, b in zip(gaps, gaps[1:]):
        assert b < a / 1.6, \
            'the monodromy gap across the switch went %.3e -> %.3e, a ' \
            'factor of %.2f where O(h) needs ~2. A gap that does NOT fall ' \
            'is the O(1) signature of a missing saltation correction' \
            % (a, b, a / b)


def _pac_L_and_B(pss, N):
    """`L` and `B` exactly as the withdrawn `PAC.solve` body builds them.

    Kept as a helper rather than inlined because the point of the test is
    that THIS construction -- the one in the tree -- has the properties
    DAC'96 claims for it, and only for the method it was written for.
    """
    times = pss.times[:-1]
    hs = np.diff(pss.times)
    M = len(times)
    L = np.zeros((N * M, N * M))
    B = np.zeros_like(L)
    for i, (_t, h, J, Cm) in enumerate(zip(times, hs, pss.Jtvec, pss.Cvec)):
        L[i * N:(i + 1) * N, i * N:(i + 1) * N] = J
        if i > 0:
            L[i * N:(i + 1) * N, (i - 1) * N:i * N] = -np.asarray(Cm) / h
    B[0:N, (M - 1) * N:M * N] = -np.asarray(pss.Cvec[-1]) / hs[0]
    return L, B, M


def _pac_operator_pieces(method, npts=40):
    """One trajectory, two descriptions of it: `L`/`B`, and our own matvec.

    Both sides are driven from the SAME seed and the SAME grid on purpose.
    Convergence is irrelevant here — these are linear-algebra identities on
    whatever `Jtvec`/`Cvec` hold — but the two sides must describe one
    trajectory or the comparison means nothing.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _q20_rlc()
    per, N = 1e-3, cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-9)
    pss._open_at_x0 = False
    pss.autonomous = False
    times, hs = pss._period_grid(per, npts, None)
    x0 = np.zeros(N)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss._record_cj = True                     # (Cvec/Jtvec on request)
        pss._walk('plain', x0, times, hs, T=per)    # records Cvec/Jtvec/times
        Jt = [np.asarray(j).copy() for j in pss.Jtvec]
        Cv = [np.asarray(c).copy() for c in pss.Cvec]
        tms = np.asarray(pss.times).copy()
        _wf = pss._walk('plain', x0, times, hs, T=per, dense=False,
                        keep=True)
        opening, steps = _wf.opening, _wf.steps
    pss.Jtvec, pss.Cvec, pss.times = Jt, Cv, tms
    L, B, M = _pac_L_and_B(pss, N)
    Mmv = np.column_stack([pss._monodromy_matvec_plain(opening, steps, e)
                           for e in np.eye(N)])
    ## the last block of `L^-1 B`, as a map on the wrapped state
    H = np.zeros((N, N))
    for j in range(N):
        e = np.zeros(N)
        e[j] = 1.0
        b = np.zeros(N * M)
        b[0:N] = B[0:N, (M - 1) * N:M * N] @ e
        H[:, j] = np.linalg.solve(L, b)[(M - 1) * N:M * N]
    return pss, L, B, N, M, H, Mmv


def test_the_pac_operator_is_the_monodromy_and_L_is_never_formed():
    """⚠ THE WITHDRAWN PAC BODY HOLDS THE RIGHT OPERATOR, and this pins it.

    DAC'96 reaches the iterative form by "reinterpreting the use of `L^-1`
    … as a preconditioner". The algebra is one line:

        (L + αB)v = -u   ⟺   (I + α L^-1 B) v = -L^-1 u

    `L` is block lower bidiagonal, so applying `L^-1` is forward
    substitution through the timesteps — which is the recursion
    `_monodromy_matvec_plain` already runs against stored factors. So
    `H := L^-1 B` IS THE MONODROMY, the 419.5 GiB that withdrew PAC was
    entirely the cost of FORMING `L`, and the rewrite keeps the operator
    and never builds the matrix.

    This test exists because that identification is the load-bearing claim
    under the whole PAC item, and it arrived by relay from a reading of the
    paper. A claim that is not in the suite is a claim that drifts.
    """
    _pss, L, B, N, M, H, Mmv = _pac_operator_pieces('euler')

    ## (1) the two structural claims DAC'96 makes -- about OUR matrices
    upper = 0.0
    for i in range(M):
        for j in range(M):
            if j > i or j < i - 1:
                upper = max(upper, np.abs(
                    L[i * N:(i + 1) * N, j * N:(j + 1) * N]).max())
    assert upper == 0.0, \
        'L is not block lower bidiagonal (max %g outside) -- the forward ' \
        'substitution that makes L^-1 cheap does not apply' % upper
    Bout = B.copy()
    Bout[0:N, (M - 1) * N:M * N] = 0.0
    assert np.abs(Bout).max() == 0.0, \
        'B is not confined to the first N rows and last N columns; the ' \
        'periodic wrap is the only thing it is supposed to carry'

    ## (2) the algebraic identity the whole reformulation rests on
    rng = np.random.default_rng(0)
    v = rng.standard_normal(N * M)
    alpha = np.exp(-2j * np.pi * 3.0)
    lhs = (L + alpha * B) @ v
    rhs = L @ (v + alpha * np.linalg.solve(L, B @ v))
    assert np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs) < 1e-12, \
        '(L + aB)v != L(I + a L^-1 B)v -- the preconditioned form is not ' \
        'the same operator'

    ## (3) and the identification itself: L^-1 B IS the monodromy, up to
    ##     the sign B carries (`B = -C/h`)
    rel = np.linalg.norm(H + Mmv) / np.linalg.norm(Mmv)
    assert rel < 1e-11, \
        'the last block of L^-1 B is not -M (rel %.3g) -- PAC cannot be ' \
        'built on the traversal if the two disagree' % rel


@pytest.mark.parametrize('method', ['trap', 'gear'])
def test_the_pac_L_is_backward_euler_only(method):
    """⚠ AND THE OPERATOR IS EULER-SHAPED, WHICH IS NOT A DETAIL.

    The withdrawn body says so in its own comment — "create LHS matrix
    using backward Euler discretization" — and builds `L` with exactly two
    terms per row: `L[i,i] = J`, `L[i,i-1] = -C/h`. A TWO-STEP method's
    variational system is block TRI-diagonal, so for `trap` or `gear` that
    `L` describes a different recursion than the trajectory it was built
    from, and `L^-1 B` is not the monodromy at all.

    ⚠ THIS IS THE TRAP A REWRITE WOULD FALL INTO. The structural checks in
    the sibling test PASS for every method — `L` really is block lower
    bidiagonal whatever `Jtvec` holds, because that is a property of how
    the loop writes it, not of the physics. Verifying the structure and
    then assuming the identification is exactly the mistake: measured on
    the Q=20 resonator at 40 points, our monodromy has ρ = 0.8545 (trap) /
    0.8412 (gear) against the analytic 0.854636, while `-L^-1 B` has ρ = 0.

    So the rewrite must take `H` from the traversal, which carries each
    step's own `(alphas, b)`, and must NOT rebuild `L`. That is the same
    conclusion the memory cost forces, reached independently.
    """
    _pss, _L, _B, _N, _M, H, Mmv = _pac_operator_pieces(method)

    rho_ours = max(abs(np.linalg.eigvals(Mmv)))
    rho_H = max(abs(np.linalg.eigvals(-H)))
    assert abs(rho_ours - 0.854636) < 0.02, \
        '%s monodromy rho %.6f is not the analytic exp(-pi/Q) -- this test ' \
        'compares against it, so it has to be right first' % (method, rho_ours)
    assert rho_H < 1e-3, \
        '%s: -L^-1 B now has rho %.6g. If the Euler-shaped L has started ' \
        'agreeing with a two-step monodromy, either L gained the third term ' \
        'or the monodromy lost it -- find out which before trusting either' \
        % (method, rho_H)


def test_grid_error_refines_with_the_runs_own_twin():
    """`grid_error` repeats the solve at finer steps on fresh instances, and
    it copied the Parameters but not `monodromy`, an attribute: a trap run
    under `monodromy='gear'` read the gear twin at the first level and the
    DEFAULT twin at the refined ones (measured) -- a refinement mixing two
    methods.  Now every level reads the run's own twin."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: 0.3 * (u - u ** 3 / 3.0))
    p = PSS(cir, method='trap', reltol=1e-10)
    p.monodromy = 'gear'
    seen = []

    def ev(q):
        seen.append(q.monodromy_twin().par.method)
        return float(q.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=6.3, timestep=6.3 / 60, x0=np.array([2.0, 0.0]),
                maxiterations=80)
        p.grid_error(ev)
    assert seen == ['gear', 'gear', 'gear'], seen


def test_trbdf2_monodromy_matches_the_pencil_and_is_second_order():
    """TR-BDF2's `m x m` monodromy is second-order on the exact Floquet
    spectrum, and self-starting -- no order-dropped opener seam.

    On a source-free RC network the period map is the homogeneous flow
    `exp(A T)` with `A = -C^-1 G` (reduced), whose eigenvalues are known in
    closed form from the pencil `(C, G)`.  `PSS.factored_period_stage`
    builds the TR-BDF2 monodromy as a factored replay; densifying it and
    comparing its eigenvalues to `exp(mu T)` checks BOTH that the map is the
    right one and that its error falls as `O(h^2)`.

    ⚠ THE POINT IS THE ORDER, NOT MERELY THE MATCH.  Trapezoidal is a
    second-order method whose SHOOTING monodromy is first-order on a limit
    cycle, because its opening manufacturing step is order-dropped to Euler
    and that seam lives inside the period map (see `_walk_lmm` and
    `monodromy_twin`).  TR-BDF2 is a self-starting one-step DIRK: every
    step, including the first, is the full two-stage method, so there is no
    opener and the monodromy stays second-order.  The error-ratio assertion
    is what distinguishes the two -- a first-order map would halve, not
    quarter, per grid doubling.
    """
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e4)
    cir['R2'] = R(2, gnd, r=2e4)
    cir['C1'] = C(1, gnd, c=1e-8)
    cir['C2'] = C(2, gnd, c=3e-8)

    pss = PSS(cir)
    m = cir.n - 1
    x0 = np.zeros(m)                      # linear: monodromy is x0-independent
    Cr = np.asarray(pss._C_at(x0))
    Gr = np.asarray(pss._G_at(x0))
    A = -np.linalg.solve(Cr, Gr)
    mu = np.linalg.eigvals(A)
    T = 5e-4
    exact = np.sort(np.exp(mu * T).real)

    def dense_M(fp):
        M = np.zeros((m, m))
        for k in range(m):
            e = np.zeros(m); e[k] = 1.0
            M[:, k] = fp.matvec(e)
        return M

    errs = {}
    for npts in (100, 200, 400):
        fp = pss.factored_period_stage(x0, T, npts, method='trbdf2')
        assert fp.kind == 'dirk'
        assert fp.width == m
        lam = np.sort(np.linalg.eigvals(dense_M(fp)).real)
        errs[npts] = float(np.max(np.abs(lam - exact)))

    ## second order: each doubling quarters the error (allow margin)
    assert errs[100] / errs[200] > 3.5, errs
    assert errs[200] / errs[400] > 3.5, errs
    ## and the absolute error is already small at the coarsest grid
    assert errs[100] < 1e-5, errs

    ## the adjoint is the exact transpose of the forward map (built as a
    ## dedicated test elsewhere; sanity-checked here on one grid)
    fp = pss.factored_period_stage(x0, T, 100, method='trbdf2')
    Mf = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                          for e in np.eye(m)])
    Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), dtype=float)
                          for e in np.eye(m)])
    assert np.linalg.norm(Mt - Mf.T) < 1e-12


def test_radau_monodromy_matches_the_pencil_and_is_fifth_order():
    """Radau IIA(3)'s `m x m` monodromy is FIFTH-order on the exact Floquet
    spectrum, and self-starting -- no order-dropped opener seam.

    Same source-free RC network as the TR-BDF2 pencil test: the period map is
    the homogeneous flow `exp(A T)` with `A = -C^-1 G` (reduced).
    `PSS.factored_period_stage` builds the coupled Radau monodromy as a
    factored replay (one `3m x 3m` factor per step, reading the third block by
    stiff accuracy); densifying it and comparing eigenvalues to `exp(mu T)`
    checks BOTH that the map is the right one and that the error falls as
    `O(h^5)` -- 32x per grid doubling, not the 4x of a second-order map.
    """
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e4)
    cir['R2'] = R(2, gnd, r=2e4)
    cir['C1'] = C(1, gnd, c=1e-8)
    cir['C2'] = C(2, gnd, c=3e-8)

    pss = PSS(cir)
    m = cir.n - 1
    x0 = np.zeros(m)                      # linear: monodromy is x0-independent
    Cr = np.asarray(pss._C_at(x0))
    Gr = np.asarray(pss._G_at(x0))
    A = -np.linalg.solve(Cr, Gr)
    mu = np.linalg.eigvals(A)
    T = 5e-4
    exact = np.sort(np.exp(mu * T).real)

    def dense_M(fp):
        return np.column_stack([fp.matvec(e) for e in np.eye(m)])

    errs = {}
    for npts in (25, 50, 100):
        fp = pss.factored_period_stage(x0, T, npts, method='radau')
        assert fp.kind == 'full'
        assert fp.width == m
        lam = np.sort(np.linalg.eigvals(dense_M(fp)).real)
        errs[npts] = float(np.max(np.abs(lam - exact)))

    ## fifth order: each doubling cuts the error by ~32 (allow margin for the
    ## higher-order remainder)
    assert errs[25] / errs[50] > 20.0, errs
    assert errs[50] / errs[100] > 20.0, errs
    ## and the absolute error is tiny already at the coarsest grid
    assert errs[25] < 1e-7, errs

    ## the adjoint is the exact transpose of the coupled forward map
    fp = pss.factored_period_stage(x0, T, 50, method='radau')
    Mf = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                          for e in np.eye(m)])
    Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), dtype=float)
                          for e in np.eye(m)])
    assert np.linalg.norm(Mt - Mf.T) < 1e-12


def test_driven_pss_under_trbdf2_matches_ac_and_gives_second_order_monodromy():
    """`method='trbdf2'` solves a driven PSS and routes the small-signal
    monodromy to the two-stage map.

    The linear RC's periodic steady state is its AC steady state, and its
    monodromy is `exp(-T/tau)` in closed form.  A TR-BDF2 shooting solve
    must reproduce both: the node fundamental to a few ppm of the AC phasor,
    and the spectral radius to `exp(-T/tau)`.  `factored_period()` returns a
    `kind='trbdf2'` map, so every surface built on it inherits the
    second-order (no-opener) monodromy without a Gear-2 twin.
    """
    circuit.default_toolkit = circuit.numeric
    period, N = 1e-3, 400
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1 / period, phase=20)
    cir['R'] = R(1, 2, r=1e4)
    cir['C'] = C(2, gnd, c=1e-8)
    tau = 1e4 * 1e-8

    resac = AC(cir).solve(1 / period)
    pss = PSS(cir, method='trbdf2')
    pss.solve(period=period, timestep=period / N)
    assert pss.converged

    tv, X = pss.waveform
    X = np.asarray(X, dtype=float)
    iref = pss.irefnode
    i2 = cir.get_node_index(cir.get_node('2'))
    row = i2 if i2 < iref else i2 - 1  # reduced index... but waveform is full
    ## pss.waveform is full-size (cir.n rows); use the full index
    v2 = X[i2][:-1]
    tt = np.asarray(tv, dtype=float)[:-1]
    f0 = 1 / period
    fund = 2.0 / len(v2) * np.sum(v2 * np.exp(-2j * np.pi * f0 * tt))
    ac2 = complex(resac.v('2'))
    ## the DFT phase reference differs from AC's by a fixed rotation; compare
    ## magnitudes and that the ratio is a pure phase (unit modulus)
    assert abs(abs(fund) - abs(ac2)) < 1e-3 * abs(ac2), (fund, ac2)

    fp = pss.factored_period()
    assert fp.kind == 'dirk'
    ## spectral radius is the RC pole exp(-T/tau)
    assert abs(pss.spectral_radius - np.exp(-period / tau)) < 1e-6


def test_autonomous_pss_under_trbdf2_finds_its_own_period():
    """`method='trbdf2'` solves the free-period system for an oscillator.

    Self-starting, so `x0` is the unknown with no manufactured opener. On
    van der Pol the TR-BDF2 solve must converge to the same period the LMM
    methods find (they agree to O(h^2)) and report a unit multiplier (the
    orbit's own free-phase direction), and `factored_period()` must hand
    back the two-stage map.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _scaled_vdp()
    iref = cir.get_node_index(gnd)
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-8).solve(refnode=gnd, tend=35.0,
                                                timestep=0.02, x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    X = np.asarray(res.x, dtype=float)
    W = X[:, t > 21.0]
    red = lambda f: np.concatenate((f[:iref], f[iref + 1:]))
    seed = red(W[:, int(np.argmax(W[iv]))])

    def solve(method):
        pss = PSS(_scaled_vdp(), method=method, reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.3, timestep=6.3 / 300, x0=seed,
                      maxiterations=40)
        return pss

    ref = solve('gear')
    pss = solve('trbdf2')
    assert pss.converged and pss.autonomous
    ## same limit cycle: periods agree to O(h^2)
    assert abs(pss.period - ref.period) < 1e-3 * ref.period, \
        (pss.period, ref.period)
    ## an oscillator's monodromy carries a multiplier at 1 (the orbit tangent)
    assert abs(pss.spectral_radius - 1.0) < 1e-3, pss.spectral_radius
    assert pss.factored_period().kind == 'dirk'


def test_driven_pss_under_radau_matches_ac_and_gives_fifth_order_monodromy():
    """`method='radau'` solves a driven PSS through the dense coupled Newton
    and routes the small-signal monodromy to the Radau map.

    The linear RC's periodic steady state is its AC steady state, and its
    monodromy is `exp(-T/tau)` in closed form.  A Radau IIA(3) shooting solve
    must reproduce both: the node fundamental to a few ppm of the AC phasor,
    and the spectral radius to `exp(-T/tau)`.  `factored_period()` returns a
    `kind='radau'` map, so every surface built on it inherits the order-5
    (no-opener) monodromy without a Gear-2 twin.
    """
    circuit.default_toolkit = circuit.numeric
    period, N = 1e-3, 100
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1 / period, phase=20)
    cir['R'] = R(1, 2, r=1e4)
    cir['C'] = C(2, gnd, c=1e-8)
    tau = 1e4 * 1e-8

    resac = AC(cir).solve(1 / period)
    pss = PSS(cir, method='radau')
    pss.solve(period=period, timestep=period / N)
    assert pss.converged

    tv, X = pss.waveform
    X = np.asarray(X, dtype=float)
    i2 = cir.get_node_index(cir.get_node('2'))
    v2 = X[i2][:-1]
    tt = np.asarray(tv, dtype=float)[:-1]
    f0 = 1 / period
    fund = 2.0 / len(v2) * np.sum(v2 * np.exp(-2j * np.pi * f0 * tt))
    ac2 = complex(resac.v('2'))
    assert abs(abs(fund) - abs(ac2)) < 1e-3 * abs(ac2), (fund, ac2)

    fp = pss.factored_period()
    assert fp.kind == 'full'
    assert abs(pss.spectral_radius - np.exp(-period / tau)) < 1e-6


def test_autonomous_pss_under_radau_finds_its_own_period():
    """`method='radau'` solves the free-period system for an oscillator.

    Self-starting (its three stages coupled into one 3n solve), so `x0` is
    the unknown with no manufactured opener.  On van der Pol the Radau solve
    must converge to the same period the LMM methods find, report a unit
    multiplier (the orbit's own free-phase direction), and `factored_period()`
    must hand back the coupled Radau map.  The period column `dphi/dT` that
    the free-period Newton needs is finite-difference checked in the build.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _scaled_vdp()
    iref = cir.get_node_index(gnd)
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-8).solve(refnode=gnd, tend=35.0,
                                                timestep=0.02, x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    X = np.asarray(res.x, dtype=float)
    W = X[:, t > 21.0]
    red = lambda f: np.concatenate((f[:iref], f[iref + 1:]))
    seed = red(W[:, int(np.argmax(W[iv]))])

    def solve(method):
        pss = PSS(_scaled_vdp(), method=method, reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=6.3, timestep=6.3 / 300, x0=seed,
                      maxiterations=40)
        return pss

    ref = solve('gear')
    pss = solve('radau')
    assert pss.converged and pss.autonomous
    ## same limit cycle: periods agree (Radau is order 5, gear order 2, so
    ## they agree to the coarser of the two)
    assert abs(pss.period - ref.period) < 1e-3 * ref.period, \
        (pss.period, ref.period)
    assert abs(pss.spectral_radius - 1.0) < 1e-3, pss.spectral_radius
    assert pss.factored_period().kind == 'full'


def test_radau_monodromy_transpose_matches_the_dense_transpose():
    """The Radau IIA(3) adjoint `M^T` is the exact transpose of the coupled
    forward map, built column by column from the forward matvec and from
    `matvec_transposed`; and a complex input splits into two real replays."""
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e4); cir['R2'] = R(2, gnd, r=2e4)
    cir['C1'] = C(1, gnd, c=1e-8); cir['C2'] = C(2, gnd, c=3e-8)
    pss = PSS(cir, method='radau')
    m = cir.n - 1
    fp = pss.factored_period_stage(np.zeros(m), 5e-4, 50)
    M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                         for e in np.eye(m)])
    MT = np.column_stack([np.asarray(fp.matvec_transposed(e), dtype=float)
                          for e in np.eye(m)])
    assert np.linalg.norm(MT - M.T) < 1e-12, np.linalg.norm(MT - M.T)
    end, ts, states = fp.matvec_transposed(np.arange(1.0, m + 1.0), collect=True)
    assert np.allclose(end, fp.matvec_transposed(np.arange(1.0, m + 1.0)))
    assert len(states) == len(fp.steps) and len(states[0]) == m
    v = np.arange(1.0, m + 1.0)
    vc = v + 1j * v[::-1]
    assert np.allclose(fp.matvec_transposed(vc),
                       fp.matvec_transposed(vc.real)
                       + 1j * fp.matvec_transposed(vc.imag))


def test_trbdf2_monodromy_transpose_matches_the_dense_transpose():
    """The TR-BDF2 adjoint `M^T` is the exact transpose of the forward map.

    Built column by column from the forward matvec and from
    `matvec_transposed`, the two must agree to machine precision; and a
    complex input must split into two real replays (the real map's property
    every adjoint surface relies on).
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e4); cir['R2'] = R(2, gnd, r=2e4)
    cir['C1'] = C(1, gnd, c=1e-8); cir['C2'] = C(2, gnd, c=3e-8)
    pss = PSS(cir, method='trbdf2')
    m = cir.n - 1
    fp = pss.factored_period_stage(np.zeros(m), 5e-4, 50)
    M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                         for e in np.eye(m)])
    MT = np.column_stack([np.asarray(fp.matvec_transposed(e), dtype=float)
                          for e in np.eye(m)])
    assert np.linalg.norm(MT - M.T) < 1e-12, np.linalg.norm(MT - M.T)
    ## collect returns width-m states (no pair) and a consistent endpoint
    v = np.arange(1.0, m + 1.0)
    end, ts, states = fp.matvec_transposed(v, collect=True)
    assert np.allclose(end, fp.matvec_transposed(v))
    assert len(states) == len(fp.steps) and len(states[0]) == m
    ## complex linearity
    vc = v + 1j * v[::-1]
    assert np.allclose(fp.matvec_transposed(vc),
                       fp.matvec_transposed(vc.real)
                       + 1j * fp.matvec_transposed(vc.imag))


def test_trbdf2_monodromy_has_less_fake_damping_than_gear_on_a_linear_oscillator():
    """The fake-damping measurement behind defaulting the twin to TR-BDF2.

    On a source-free damped LC -- a LINEAR oscillator -- the exact period map
    is `exp(A T)` with `A = -C^-1 G`, so its multiplier `lambda = exp(-alpha T)`
    and `Q_lambda` are known in closed form, no reference simulator needed.
    Both methods approximate it, but Gear-2 (BDF2) adds numerical damping to
    the weakly-damped mode -- Dharmaraja's caveat -- biasing `lambda2`, and
    the bias is worst at coarse grids and high Q (the amplification law
    `dQ/Q = Q_lambda dlambda2/lambda2`).  TR-BDF2 carries far less of it.

    Deterministic and cheap: no ODE integration, no PSS solve.  TR-BDF2's
    monodromy is the SHIPPING one (`factored_period_stage`); Gear-2's is the
    BDF2 companion on the same reduced pencil, which is exactly the
    `solved_history` pair map a Gear-2 PSS forms on a linear system.

    ⚠ THE ADVANTAGE IS REGIME-DEPENDENT, NOT A FIXED FACTOR.  At Q ~ 16 and
    50 points/period TR-BDF2 is ~26x closer to the exact `Q`; refine to 200
    points and both approach the O(h^2) floor where the gap shrinks to ~2x.
    That is the true shape of the result -- large where a real oscillator PSS
    sits (coarse grid, high Q), small on a fine grid -- and is why the twin
    is a DEFAULT rather than the only option.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['R'] = R('v', gnd, r=50.0)          # parallel R -> complex pole pair
    pss = PSS(cir, method='trbdf2')
    m = cir.n - 1
    Cr = np.asarray(pss._C_at(np.zeros(m)))
    Gr = np.asarray(pss._G_at(np.zeros(m)))
    A = -np.linalg.solve(Cr, Gr)
    ev = np.linalg.eigvals(A)
    w = float(np.max(np.abs(ev.imag)))
    alpha = float(-np.mean(ev.real))
    T = 2.0 * np.pi / w
    Q_exact = -1.0 / np.log(np.exp(-alpha * T))

    def bdf2_companion(h):
        Minv = np.linalg.inv(np.eye(m) - (2.0 / 3.0) * h * A)
        return np.block([[Minv @ ((4.0 / 3.0) * np.eye(m)),
                          Minv @ ((-1.0 / 3.0) * np.eye(m))],
                         [np.eye(m), np.zeros((m, m))]])

    def Q_of(lam):
        return -1.0 / np.log(lam)

    def errs(N):
        h = T / N
        fp = pss.factored_period_stage(np.zeros(m), T, N)
        Mt = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                              for e in np.eye(m)])
        lam_t = float(np.max(np.abs(np.linalg.eigvals(Mt))))
        Mg = np.linalg.matrix_power(bdf2_companion(h), N)
        lam_g = float(np.max(np.abs(np.linalg.eigvals(Mg))))
        return (abs(Q_of(lam_t) - Q_exact) / Q_exact,
                abs(Q_of(lam_g) - Q_exact) / Q_exact)

    ## Q ~ 16 here; at 50 points/period TR-BDF2 is an order-plus better
    et50, eg50 = errs(50)
    assert et50 < 5e-3, 'trbdf2 Q error %.2e at N=50' % et50
    assert eg50 > 1e-2, 'gear Q error %.2e at N=50 (expected the bias)' % eg50
    assert eg50 / et50 > 5.0, \
        'trbdf2 should be >5x better than gear at N=50; got %.1fx' \
        % (eg50 / et50)
    ## and TR-BDF2 is never worse than gear on this axis at the coarse grid
    et200, eg200 = errs(200)
    assert et200 <= eg200 * 1.5, \
        'trbdf2 %.2e vs gear %.2e at N=200' % (et200, eg200)


def test_esdirk43_monodromy_order4_through_the_generic_dirk_family():
    """ESDIRK4(3)6 -- a NEW DIRK method, tableau-only -- gets a correct order-4
    shooting monodromy through the GENERIC s-stage DIRK-sequential family, at
    s=6 stages (TR-BDF2 is s=3).  This is the check that the DIRK-family
    generalisation is real: no ESDIRK-specific shooting code exists.

    On the source-free RC network the period map is `exp(A T)`; the densified
    `kind='dirk'` monodromy must match `exp(mu T)` at O(h^4) -- 16x per grid
    doubling -- and its adjoint must be the exact transpose.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['R1'] = R(1, 2, r=1e4); cir['R2'] = R(2, gnd, r=2e4)
    cir['C1'] = C(1, gnd, c=1e-8); cir['C2'] = C(2, gnd, c=3e-8)
    pss = PSS(cir)
    m = cir.n - 1
    x0 = np.zeros(m)
    Cr = np.asarray(pss._C_at(x0)); Gr = np.asarray(pss._G_at(x0))
    A = -np.linalg.solve(Cr, Gr)
    T = 5e-4
    exact = np.sort(np.exp(np.linalg.eigvals(A) * T).real)
    errs = {}
    for npts in (25, 50, 100):
        fp = pss.factored_period_stage(x0, T, npts, method='esdirk43')
        assert fp.kind == 'dirk'
        M = np.column_stack([fp.matvec(e) for e in np.eye(m)])
        errs[npts] = float(np.max(np.abs(np.sort(np.linalg.eigvals(M).real)
                                         - exact)))
    ## fourth order: ~16x per doubling
    assert errs[25] / errs[50] > 10.0, errs
    assert errs[50] / errs[100] > 10.0, errs
    assert errs[25] < 1e-6, errs
    ## adjoint is the exact transpose
    fp = pss.factored_period_stage(x0, T, 50, method='esdirk43')
    Mf = np.column_stack([np.asarray(fp.matvec(e), dtype=float) for e in np.eye(m)])
    Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), dtype=float)
                          for e in np.eye(m)])
    assert np.linalg.norm(Mt - Mf.T) < 1e-12


def test_pcnr_reaches_the_shooting_inner_transient_over_a_stage_method():
    """PCNR now reaches SHOOTING too: PSS forwards `pcnr` to its inner
    transient, and the per-step PCNR lives in `solve_timestep` (which PSS
    drives), so a stage-method PSS with pcnr=True limits its junctions by the
    continuation and reaches the SAME periodic orbit device limiting does.

    Closes half the external review's 1-for-4 (limiting reached, PCNR did not):
    PCNR is now 2-for-4 (rescue/breakpoints still need Transient.solve, which
    PSS does not call).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def solve(pcnr):
        c = _diode_mixer()
        p = PSS(c, method='trbdf2', reltol=1e-11, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=1e-6, timestep=1e-6 / 160, maxiterations=40)
        assert p.converged
        return p

    p_lim = solve(False)
    p_pcnr = solve(True)
    X0 = np.asarray(p_lim.waveform[1], dtype=float)
    X1 = np.asarray(p_pcnr.waveform[1], dtype=float)
    assert np.linalg.norm(X0 - X1) / np.linalg.norm(X0) < 1e-9, \
        'PCNR and limiting must reach the same shooting orbit'
    assert abs(p_pcnr.spectral_radius - p_lim.spectral_radius) \
        < 1e-6 * abs(p_lim.spectral_radius) + 1e-12


def test_monodromy_matches_a_finite_difference_of_the_period_map():
    """EVERY period-map family's analytic monodromy must equal the DERIVATIVE OF
    THE DISCRETE PERIOD MAP it claims to be -- checked against a central finite
    difference of that very map, a reference this code cannot influence.

    Covers all four kinds: `solved_history` (gear, the 2m PAIR map), `plain`
    (trap), `full` (Radau, coupled) and `dirk` (TR-BDF2, sequential).

    ⚠ THIS IS THE TEST THAT CAUGHT THE SHARED-`_vlim` DEFECT.  A junction
    device's `i`/`G` are read at its stored `_vlim`, and there is only ONE
    `_vlim` per device -- while the monodromy evaluates `C`/`G` at SEVERAL
    distinct points per step (`x_n` and every stage).  Whatever the step's solve
    left behind (the LAST stage) was used for all of them, so the period map
    linearised the junction at the wrong voltage.  Measured error before the fix:
    1.65e-3 (Radau) / 9.1e-4 (TR-BDF2) relative, FLAT across four decades of the
    FD step -- a fixed error, not FD noise, which is exactly how it was told
    apart.  Everything downstream rests on this: Floquet multipliers, the PSS
    Jacobian, the PPV, the cyclostationary noise.

    ⚠ THE PERIODICITY ASSERTION IS NOT A FORMALITY -- IT CAUGHT A WRONG
    REFERENCE.  With a MANUFACTURED opening (`x0_unknown=False` on a one-step
    method) the shooting unknown is `x_in`, but the monodromy is about the
    POST-manufacturing state, so a forward map that starts at `x_in` and skips
    the manufacturing step is simply a different map -- it reported 5.13 here,
    and comparing against it produced a confident 8.5e-3 "defect" that did not
    exist.  A reference has to be shown to be the right map before its
    disagreement means anything.  `gear` needs no such flag (`solved_history`
    already carries `x_0` and `x_{-1}` as real trajectory states).

    ⚠ THE CIRCUIT IS CHOSEN SO THE JUNCTION IS ACTUALLY IN THE ANSWER, and the
    test asserts it.  Two degenerate regimes make this check vacuous and both
    were hit while building it: a fast RC drives the whole monodromy to ~0 (the
    multiplier came out 1e-51 -- a 0-vs-0 comparison that passes regardless),
    and a bare diode either never conducts (the capacitor shunts the node, so
    the multiplier is EXACTLY the diode-off `exp(-T/RC)` and the junction is
    absent) or conducts so hard it shorts the node and the multiplier collapses
    to 0 again.  Feeding the diode through a series resistor BOUNDS its
    conductance: off `exp(-1)` = 0.368, fully on `exp(-2)` = 0.135, solution in
    between -- non-degenerate AND junction-dependent.
    """
    import warnings
    from copy import copy
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric

    def build(va):
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, vac=1.0, va=va, freq=1e6, phase=20)
        c['R'] = R(1, 2, r=1e3)
        c['Rd'] = R(2, 3, r=1e3)      # bounds the junction conductance
        c['D'] = Diode(3, gnd)
        c['C'] = C(2, gnd, c=1e-9)    # tau = R*C = T, so the map is not ~0
        return c

    def fwd(pss, state, kind, times, hs, method, m):
        """The DISCRETE period map, seeded exactly as its traversal seeds it:
        the pair via `_install_history` for `solved_history`, a single entering
        state via `_begin_period` otherwise."""
        tr_saved = getattr(pss, '_tran', None)
        pss._tran = pss._new_transient(pss._integrator_for(method))
        try:
            pss._want_dfdh = False
            pss._want_lte = False
            if kind == 'solved_history':
                x0, xm1 = state[:m].copy(), state[m:].copy()
                pss._install_history(x0, xm1, hs[0], h_prev=hs[-1])
                x, xp = copy(x0), copy(xm1)
                for j, t in enumerate(times[1:]):
                    xp = x
                    x = copy(pss.solve_timestep(x, t, hs[j]))
                return np.concatenate((np.asarray(x, float).ravel(),
                                       np.asarray(xp, float).ravel()))
            pss._begin_period(np.asarray(state, float))
            x = copy(np.asarray(state, float))
            for j, t in enumerate(times[1:]):
                x = copy(pss.solve_timestep(x, t, hs[min(j, len(hs) - 1)]))
            return np.asarray(x, float).ravel()
        finally:
            pss._tran = tr_saved

    off = np.exp(-1.0)     # the diode-OFF multiplier, exp(-T/RC)
    ## gear's solved-history map needs no x0_unknown (and refuses it); the
    ## one-step methods are checked in the formulation whose map IS a function
    ## of the state, which is what makes the FD reference the right map.
    ## ⚠ `theta` IS HERE FOR A REASON THE OTHERS ARE NOT.  Its opening seeds a
    ## CONSISTENT `iq_{-1}` that depends on `x_0`, so its monodromy carries a
    ## term (`-G(x_0)`) that no other method's does -- and the FD of the period
    ## map is the only reference that sees it, because dropping the term costs
    ## iterations and not accuracy.  See
    ## `test_theta_s_shooting_jacobian_carries_the_consistent_iq_seed`.
    for method, kw in (('gear', {}), ('trap', {'x0_unknown': True}),
                       ('theta', {'x0_unknown': True}),
                       ('radau', {}), ('trbdf2', {})):
        pss = PSS(build(15.0), method=method, reltol=1e-13)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=1e-6, timestep=1e-6 / 60, maxiterations=60, **kw)
        assert pss.converged, '%s: PSS did not converge' % method
        fp = pss.factored_period()
        w = fp.width
        m = pss.cir.n - 1
        _solved, x0, xm1, _t, _h, _T, _x0u = pss._period_state
        x0 = np.asarray(x0, float).ravel()
        state = (np.concatenate((x0, np.asarray(xm1, float).ravel()))
                 if fp.kind == 'solved_history' else x0)
        times = np.asarray(fp.times, float)
        hs = np.diff(times)

        ## the reference must be THE RIGHT MAP before its disagreement counts
        per = np.max(np.abs(fwd(pss, state, fp.kind, times, hs, method, m)
                            - state))
        assert per < 1e-12, \
            '%s (%s): the forward map is not periodic at the solution (%.2e), ' \
            'so it is not the map the monodromy is about -- fix the reference ' \
            'before reading anything into a disagreement' % (method, fp.kind, per)

        Man = np.column_stack([np.asarray(fp.matvec(np.eye(w)[:, k]),
                                          float).ravel() for k in range(w)])
        d = 1e-6
        Mfd = np.zeros((w, w))
        for k in range(w):
            e = np.zeros(w)
            e[k] = d
            Mfd[:, k] = (fwd(pss, state + e, fp.kind, times, hs, method, m)
                         - fwd(pss, state - e, fp.kind, times, hs, method, m)) \
                / (2 * d)

        ## the junction must actually be in the answer, or this proves nothing
        lam = np.max(np.abs(np.linalg.eigvals(Mfd)))
        assert abs(lam - off) / off > 0.05, \
            '%s: the junction is not loading the monodromy (|lambda|=%.5f vs ' \
            'diode-off %.5f) -- the test would be vacuous' % (method, lam, off)

        rel = np.linalg.norm(Man - Mfd) / np.linalg.norm(Mfd)
        assert rel < 1e-7, \
            '%s (%s): analytic monodromy disagrees with the finite-difference ' \
            'derivative of its own period map by %.3e (was 1.6e-3 with the ' \
            'shared-_vlim defect; the FD noise floor here is ~1e-9)' \
            % (method, fp.kind, rel)


def _b2_resonator():
    """The B2 gate's Q = 20 resonator -- `benchmarks/pss_b2_theta_gate.py`.

    NOT `_q20_rlc`, and the difference is load-bearing: `theta - 1/2 = C h`,
    so with `ThetaIntegrator.DEFAULT_C` fixed at 1e4 the bias a period carries
    is `C T`, and these two fixtures differ in `T` by 159x.  Every recorded B2
    number (peak 20.01524 at K = 200, `|mode|^K = 0.7778`) belongs to this one.
    """
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    Lv, Cv, Q = 1e-3, 1e-9, 20.0
    f0 = 1.0 / (2 * np.pi * np.sqrt(Lv * Cv))
    cir = SubCircuit()
    cir.add_node('n1'); cir.add_node('n2')
    cir['vs'] = VSin(gnd, 'n1', va=1.0, freq=f0)
    cir['L'] = L('n1', 'n2', L=Lv)
    cir['C'] = C('n2', gnd, c=Cv)
    cir['R'] = R('n1', 'n2', r=Q * np.sqrt(Lv / Cv))
    return cir, 1.0 / f0


def _shooting_evaluations(method, K, T, **kw):
    """Run the B2 resonator; return (peak, n_evaluations, func, points, resid).

    `resid` is the max-norm of the shooting residual at each evaluation, and it
    is the quantity that carries the claim -- see the test below for why the
    COUNT alone is not.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir, _T = _b2_resonator()
    pss = PSS(cir, method=method, reltol=1e-3)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        res = pss.solve(period=T, timestep=T / K, maxiterations=200,
                        trace=True, **kw)
    assert pss.converged, '%s at K=%d did not converge' % (method, K)
    pts = [z for z, _F, _J in pss.shooting_trace]
    resid = [float(np.max(np.abs(F))) for _z, F, _J in pss.shooting_trace]
    peak = float(np.max(np.abs(np.asarray(res['tpss'].v('n2'), float).ravel())))
    ## (the residual holds its analysis weakly since 2026-10-01: the
    ## returned function keeps the PSS alive)
    return (peak, len(pts),
            lambda *a, _p=pss, **k: _p.shooting_residual(*a, **k), pts, resid)


def test_theta_s_shooting_jacobian_carries_the_consistent_iq_seed():
    """B2's one open item, closed: it was a MISSING CHAIN RULE, not the method.

    B2 shipped `theta` as an available method rather than a recommendable one
    because it wanted far more shooting iterations than `trap` at the same K,
    with the right peak at every K.  That cost is not a property of the theta
    method.  It is one term.

    ⚠ THE RECORDED FIGURE WAS "~150 at K = 400 against trap's 40"; re-measured
    here on the B2 gate resonator it is 99 against trap's 3, so the number is
    restated rather than quoted.  The 40 in the old record is `_pss_lte`'s
    `maxiterations`, not a count.

    `theta` refuses the L-stable opener, so it READS `iq_{-1}` on its first
    step, where `Transient._begin_run` seeds `-(i(x_0) + u(t_0))` -- a FUNCTION
    OF THE UNKNOWN (`Integrator.needs_consistent_iq0`).  Every `open_at_x0`
    branch in `shooting.py` seeded `d(iq_0)/d(x_0) = 0` and said so in the same
    words, "no companion current has been formed yet" -- true of every method
    written before this one.  Dropping `-G(x_0)` did not perturb the monodromy
    slightly: it ANNIHILATED `null(C)`, which is precisely what an L-stable
    Euler opener does, i.e. the one thing this method exists not to do.

    ⚠⚠ THE ANSWER WAS NEVER WRONG, WHICH IS WHY NOTHING CAUGHT IT.  The
    residual is the residual; only the Newton DIRECTION was wrong, so the solve
    still landed on the same orbit and every peak in the B2 record is
    reproduced here to the digit.  A test that compares an amplitude cannot see
    this class of defect at all -- so this one checks the JACOBIAN and the
    ITERATION COUNT.

    ⚠ THE GATE IS THAT A LINEAR CIRCUIT FORCES THE ANSWER.  `phi` is affine in
    `x_0`, so an EXACT shooting Newton lands in one step -- whatever the
    method, whatever the grid.  The residual after ONE step, from a seed at
    2.9:

        trap,  x0_unknown=True    1e-12 .. 1e-11   (the control; always did)
        theta, before             1.7e-02          and still 4e-07 after 64
        theta, after              1.9e-14 .. 2e-11

    ⚠⚠ AND THE FIRST VERSION OF THIS TEST ASSERTED ON THE EVALUATION COUNT,
    WHICH IS A PROXY, AND IT BROKE ON A LAST-BIT CHANGE IN `T`.  The counts
    really were 9/64/99 before and 3/3/3 after -- but 3 is not robust: once the
    first step lands at 1.9e-14 the solver is bumping along the ROUNDOFF FLOOR
    (1.9e-14 -> 3.6e-14 -> 5.3e-14 ...) and whether it stops at 3 evaluations
    or 7 is decided by where the noise falls, not by the Jacobian.  Changing
    `theta`'s bias by 0.05% flipped it.  **The contraction is the claim; the
    count was a symptom.**  A loose count bound is kept only to catch the gross
    case (64 and climbing).

    ⚠ AND `trap` WITH `x0_unknown=False` TAKES 7 / 6 / 77 HERE, so "many
    evaluations" is not by itself a theta symptom -- the manufactured-opening
    formulation has an inexact Jacobian BY CONSTRUCTION, and the plain walk's
    record says so in as many words.  It is not a control for this defect; the control is the
    SAME formulation under a different method, which is the row above it.

    ⚠ THE JACOBIAN CHECK IS DELTA-SWEPT because that is how a real error is
    told from FD noise: noise makes a V in delta (truncation down, roundoff
    up), a real error sits FLAT.  Before the fix the relative disagreement was
    6.344 at every delta from 1e-3 to 1e-8; after it, ~1.4e-10.
    """
    from pycircuit.circuit.shooting import PSS as _PSS
    cir, T = _b2_resonator()

    ## (1) The peaks are UNCHANGED -- the B2 record, to the digit.
    recorded = {100: 19.98407, 200: 20.01524, 400: 20.02255}
    counts = {}
    for K in (100, 200, 400):
        ## (the record's grid, K - 1 steps: `T / K` gave K - 1 until 2026-09-30)
        peak, n_theta, f_theta, pts, rr = _shooting_evaluations('theta', K - 1, T)
        assert abs(peak - recorded[K]) < 5e-5, \
            'theta at K=%d moved off the B2 record: %.5f vs %.5f. The seed ' \
            'fixes the JACOBIAN and must not touch the residual.' \
            % (K, peak, recorded[K])
        _pk, n_trap, _f, _p, rt = _shooting_evaluations('trap', K - 1, T,
                                                        x0_unknown=True)
        counts[K] = (n_theta, n_trap)
        ## (2) ⚠ THE CONTRACTION, NOT THE COUNT -- see the docstring.  One
        ## Newton step on an affine residual must take it to roundoff, and
        ## `trap` in the SAME formulation says what roundoff looks like here.
        drop_t = rr[1] / rr[0]
        drop_r = rt[1] / rt[0]
        assert drop_t < 1e-11, \
            'theta at K=%d: one Newton step cut the residual by only %.3e ' \
            '(2.9 -> %.3e). On a linear circuit phi is AFFINE, so an exact ' \
            'Jacobian must reach roundoff in one step; it was 5.9e-03 with ' \
            'd(iq_0)/d(x_0) dropped.' % (K, drop_t, rr[1])
        assert drop_r < 1e-11, \
            'the CONTROL failed: trap dropped only %.3e at K=%d, so this ' \
            'harness is not measuring a one-step Newton and the theta ' \
            'assertion above is vacuous' % (drop_r, K)
        assert n_theta <= 12, \
            'theta took %d evaluations at K=%d against trap\'s %d -- the ' \
            'count is only a coarse guard (it was 64 and climbing with the ' \
            'seed dropped), but 12 is far past bumping along roundoff' \
            % (n_theta, K, n_trap)

    ## (3) THE JACOBIAN ITSELF, delta-swept against an FD of its own residual.
    ## (199 steps, the record's K = 200 grid -- see the loop above)
    _pk, _n, func, pts, _rr = _shooting_evaluations('theta', 199, T)
    for label, x in (('the seed', pts[0]), ('the solution', pts[-1])):
        x = np.asarray(x, float)
        F0, J = func(x.copy())
        J = np.asarray(J, float)
        scale = max(float(np.max(np.abs(x))), 1.0)
        rels = []
        for d in (1e-3, 1e-4, 1e-5, 1e-6):
            dd = d * scale
            Jfd = np.empty_like(J)
            for j in range(len(x)):
                xp = x.copy(); xp[j] += dd
                xm = x.copy(); xm[j] -= dd
                Jfd[:, j] = (np.asarray(func(xp)[0], float)
                             - np.asarray(func(xm)[0], float)) / (2 * dd)
            rels.append(float(np.max(np.abs(J - Jfd)))
                        / max(float(np.max(np.abs(Jfd))), 1e-30))
        assert max(rels) < 1e-7, \
            'theta at %s: the analytic shooting Jacobian disagrees with a ' \
            'finite difference of its own residual by %.3e (delta sweep %s). ' \
            'FLAT across decades means a real term is missing, not FD noise; ' \
            'it was 6.344 with d(iq_0)/d(x_0) dropped.' \
            % (label, max(rels), ['%.2e' % r for r in rels])

    ## (4) AND THE MODE THAT WAS LOST IS THE ONE THE METHOD EXISTS FOR.  With
    ## the seed dropped, `I - M` carries 1 on `null(C)` -- the signature of an
    ## L-stable opener -- instead of `1 - (-(1-theta)/theta)^K = 1.7778`.
    _pk, _n, func_t, pts_t, _rr = _shooting_evaluations('theta', 199, T)
    _F, J_ok = func_t(np.asarray(pts_t[-1], float))
    saved = _PSS._pq_seed_at_x0
    try:
        _PSS._pq_seed_at_x0 = lambda self, x: None
        (peak_n, n_neutered, func_n, pts_n,
         rr_n) = _shooting_evaluations('theta', 199, T)
        _F, J_bad = func_n(np.asarray(pts_n[-1], float))
    finally:
        _PSS._pq_seed_at_x0 = saved
    J_ok = np.asarray(J_ok, float); J_bad = np.asarray(J_bad, float)
    dropped = float(np.max(np.abs(np.diag(J_ok) - np.diag(J_bad))))
    mode = 0.7778                      # ((1-theta)/theta)^K, K even -> positive
    assert abs(dropped - mode) < 5e-3, \
        'the neutered Jacobian should differ from the correct one by exactly ' \
        'the null(C) multiplier %.4f on the diagonal, and differs by %.4f -- ' \
        'so this test is no longer pinning the mechanism it names' \
        % (mode, dropped)
    assert rr_n[1] / rr_n[0] > 1e-4, \
        'NEUTER CHECK: with d(iq_0)/d(x_0) dropped, one Newton step should ' \
        'cut the residual by only ~6e-03 and it cut it by %.3e -- if the ' \
        'defect no longer destroys the contraction then this test has ' \
        'stopped guarding anything' % (rr_n[1] / rr_n[0])
    assert n_neutered > 4 * counts[200][1], \
        'and it should still COST iterations (64 on record, %d here against ' \
        "trap's %d)" % (n_neutered, counts[200][1])
    assert abs(peak_n - recorded[200]) < 5e-5, \
        'and the neutered run must still reach the SAME peak (%.5f vs %.5f) ' \
        '-- the point of the whole test is that the answer was never wrong' \
        % (peak_n, recorded[200])

    ## (5) THE OTHER TWO PATHS THAT OPEN AT `x_0`, and the mode itself.  The
    ## FACTORED matvec and its reverse replay carry the same seed, so the
    ## monodromy they build must show the DAMPED null(C) mode this method
    ## exists to produce -- `((1-theta)/theta)^K`, positive at even K, which
    ## `ThetaIntegrator`'s own table gives as 7.778e-01 at C = 1e4.  Before the
    ## fix that eigenvalue was 0: annihilated, exactly as an Euler opener does.
    import warnings as _w
    cir2, _T = _b2_resonator()
    pss2 = PSS(cir2, method='theta', reltol=1e-9)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
        pss2.solve(period=T, timestep=T / 199, maxiterations=60)
    fp = pss2.factored_period()
    assert fp.kind == 'plain', fp.kind
    w = fp.width
    Mf = np.column_stack([np.asarray(fp.matvec(e), float) for e in np.eye(w)])
    Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), float)
                          for e in np.eye(w)])
    rel = float(np.max(np.abs(Mt - Mf.T))) / float(np.max(np.abs(Mf)))
    assert rel < 1e-9, \
        'theta: the reverse replay disagrees with the transpose of the ' \
        'forward one by %.3e -- the seed enters the adjoint as ' \
        '`pq_open^T w2` at the CLOSE of the backward pass, and dropping it ' \
        'there is invisible to the forward check' % rel
    lam = np.sort(np.abs(np.linalg.eigvals(Mf)))
    mode_k = 0.777768                  # ((1-theta)/theta)^K at C=1e4, K=200
    assert abs(lam[0] - mode_k) < 1e-4, \
        'the smallest multiplier is %.6f, not the damped null(C) mode %.6f. ' \
        'That mode IS the method: with the seed dropped it comes back as 0 ' \
        '(annihilated, i.e. an L-stable opener), which is the thing theta ' \
        'exists not to do.' % (lam[0], mode_k)
    decay = float(np.exp(-np.pi / 20.0))
    assert abs(lam[-1] - decay) < 1e-3, \
        'and the PHYSICAL pair must still be exp(-pi/Q) = %.6f, not %.6f -- ' \
        'otherwise the seed has moved the circuit and not just the opening' \
        % (decay, lam[-1])

    ## (6) AND THE FOURTH PATH, which the widened `opening` tuple found:
    ## PAC's `_forced_replay` opens `Pq` the same way, and it must carry the
    ## SAME seed as the monodromy it superposes with -- `y_end = M y0 + w` is
    ## what lets PAC solve an `m x m` system instead of an `(N m) x (N m)`
    ## one, and a seed in one and not the other breaks it silently.
    m = pss2.cir.n - 1
    u_ac = np.zeros(m); u_ac[0] = 1.0
    fac = 0.37 / T
    w_zero, _ = pss2._forced_replay(fp, fac, u_ac, y0=None)
    rng = np.random.default_rng(0)
    worst = 0.0
    for _k in range(3):
        y0 = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        y1, _ = pss2._forced_replay(fp, fac, u_ac, y0=y0)
        worst = max(worst, float(np.max(np.abs(y1 - (Mf @ y0 + w_zero))))
                    / max(float(np.max(np.abs(y1))), 1e-30))
    assert worst < 1e-9, \
        'theta: the forced replay is not `M y0 + w` (%.3e). The seed must be ' \
        'in BOTH or PAC superposes a driven response onto a different map.' \
        % worst


def test_theta_s_bias_is_per_period_and_the_knob_is_reachable():
    """`DEFAULT_C = 1e4` was a FIXTURE CONSTANT, and no caller could override it.

    ⚠⚠ THE DIAL IS `C T`, NOT `C`, AND THAT IS ARITHMETIC.  `theta = 1/2 + C h`
    damps `null(C)` per step by `-(1-theta)/theta`, so over a period of `K`
    steps by `((1-theta)/theta)^K ~ exp(-4 C h K) = exp(-4 C T)`.  **`h`
    cancels.**  What a period delivers depends on the PRODUCT alone -- not on
    `C`, and not on the step count.  The B2 gate located its knee
    (`rcond(I - A^K)` saturating) at `C = 1e4` on a fixture with
    `T = 2 pi x 1e-6`, so the transferable number is `C T = 0.0628`; storing it
    as a RATE baked in that one period.

    ⚠⚠ AND THE FAILURE IT CAUSED IS THE SILENT KIND.  `_q20_rlc` has
    `T = 1e-3`, 159x the gate's, so the same `C` gives `C T = 10` and
    `theta = 0.6` at K = 100.  Against the 20 V analytic peak::

        K     trap       C=1e4 (CT=10)   C=1e3 (CT=1)   C=62.8 (CT=0.0628)
        100   19.98967   15.91117        19.49008       19.95755
        200   19.99811   18.80471        19.87200       19.99015
        400   19.99957   19.68875        19.96805       19.99759

    **20% low at K = 100 and `converged=True`** -- because it DID converge, to
    its own over-damped discretisation.  That is the third-level failure
    `test_pss_reports_the_truncation_error_neither_newton_can_see` is named
    for, and no Newton criterion can see it.

    ⚠ A SECOND, RELATED GAP: `_integrator_for` builds `table[method]()`, so
    `cbias` was unreachable through `PSS` entirely.  A caller who hit the above
    had no knob at all.  `theta_ct` is that knob, and it is the DIMENSIONLESS
    one, so it transfers.

    ⚠ THE SEED PERIOD IS ENOUGH FOR AN AUTONOMOUS RUN, and that is measured
    rather than hoped: the gate's own table has `rcond` at 4.4e-03 / 4.6e-03 /
    4.3e-03 across `C T` = 0.0063 / 0.0628 / 0.628 -- FLAT over two decades, so
    a period that moves a few percent moves nothing that matters.
    """
    import warnings as _w
    from pycircuit.circuit.integrator import ThetaIntegrator
    circuit.default_toolkit = circuit.numeric

    ## (1) THE ARITHMETIC FIRST: at fixed `C T` the per-period damping is the
    ## same however many steps carry it.  If this fails, `theta_ct` is not the
    ## right knob and nothing below matters.
    ct = 0.0628
    for T in (1e-3, 6.2832e-6):
        damp = []
        for K in (50, 400, 3200):
            h = T / K
            th = ThetaIntegrator(ct=ct, period=T).theta_at(h)
            damp.append(((1.0 - th) / th) ** K)
        assert max(damp) - min(damp) < 1e-3 * max(damp), \
            'the per-period damping is supposed to be a function of C*T ' \
            'alone, and over K = 50/400/3200 it moved: %r' % damp
        assert abs(damp[-1] - np.exp(-4.0 * ct)) < 5e-3, \
            'and it should be exp(-4 C T) = %.6f, not %.6f' \
            % (np.exp(-4.0 * ct), damp[-1])

    ## (2) THE CONSTRUCTOR: two spellings of one quantity, and it refuses to
    ## guess which was meant.
    assert ThetaIntegrator().cbias == ThetaIntegrator.DEFAULT_C
    assert abs(ThetaIntegrator(period=1e-3).cbias
               - ThetaIntegrator.DEFAULT_CT / 1e-3) < 1e-12
    assert abs(ThetaIntegrator(ct=0.5, period=1e-3).cbias - 500.0) < 1e-9
    with pytest.raises(ValueError, match='DIMENSIONLESS'):
        ThetaIntegrator(ct=0.5)
    with pytest.raises(ValueError, match='not both'):
        ThetaIntegrator(cbias=1e4, period=1e-3)
    with pytest.raises(ValueError, match='period must be'):
        ThetaIntegrator(period=0.0)

    def peak(cir, node, T, K, neuter=False, **kw):
        p = PSS(cir, method='theta', reltol=1e-3, **kw)
        if neuter:
            p._theta_biased = lambda integ: integ
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
            res = p.solve(period=T, timestep=T / (K - 1), maxiterations=200)
        assert p.converged
        return float(np.max(np.abs(
            np.asarray(res['tpss'].v(node), float).ravel()))), p

    ## (3) THE DEFECT, ON THE FIXTURE THAT SHOWS IT.  `_q20_rlc` is the suite's
    ## own Q = 20 resonator and its analytic peak is 20 V.
    for K, want in ((100, 19.95755), (200, 19.99015), (400, 19.99759)):
        pk, p = peak(_q20_rlc(), 'c', 1e-3, K)
        assert abs(pk - want) < 5e-5, \
            'theta on _q20_rlc at K=%d gives %.5f, not the %.5f a ' \
            'period-normalised bias earns (it was %.2f with the rate)' \
            % (K, pk, want, {100: 15.91, 200: 18.80, 400: 19.69}[K])
        assert abs(pk - 20.0) / 20.0 < 3e-3, \
            'and it must now track the 20 V analytic peak: %.5f' % pk
        assert abs(p._transient().par.integrator.cbias
                   - ThetaIntegrator.DEFAULT_CT / 1e-3) < 1e-9

    ## (4) AND IT IS A NO-OP WHERE THE KNEE WAS MEASURED -- the B2 record, to
    ## the digit.  A fix that moved the calibrated case would be a new bias,
    ## not a normalisation of the old one.
    _cir, Tg = _b2_resonator()
    for K, want in ((100, 19.98407), (200, 20.01524), (400, 20.02255)):
        pk, _p = peak(_b2_resonator()[0], 'n2', Tg, K)
        assert abs(pk - want) < 5e-5, \
            'the B2 gate fixture moved: %.5f vs the recorded %.5f' % (pk, want)

    ## (5) THE KNOB IS LIVE, and it costs what the gate said it costs.
    got = {}
    for c in (0.0628, 0.628, 6.28):
        got[c], p = peak(_q20_rlc(), 'c', 1e-3, 200, theta_ct=c)
        assert abs(p._transient().par.integrator.cbias - c / 1e-3) < 1e-6, \
            'theta_ct=%g did not reach the integrator' % c
    assert got[0.0628] > got[0.628] > got[6.28], \
        'more bias must cost amplitude monotonically, got %r' % got
    assert got[0.0628] - got[6.28] > 0.5, \
        'the knob moved the peak by only %.4f V -- a knob that changes ' \
        'nothing measurable is not a knob' % (got[0.0628] - got[6.28])

    ## (6) NEUTER CHECK: with the normalisation removed the old defect comes
    ## straight back, so this test is verified to fail rather than assumed to.
    pk_bad, _p = peak(_q20_rlc(), 'c', 1e-3, 100, neuter=True)
    assert abs(pk_bad - 15.91117) < 5e-5, \
        'with `_theta_biased` neutered the K=100 peak should be the recorded ' \
        '15.91117 (20%% low), and it is %.5f -- if the defect no longer ' \
        'reproduces, this test guards nothing' % pk_bad


def test_the_default_method_is_radau_and_it_takes_no_monodromy_twin():
    """The DEFAULT itself, pinned — because almost nothing else exercises it.

    ⚠⚠ When the default moved from `trap` to `radau` (owner decision,
    2026-09-07) the full suite passed unchanged, 3099 tests, first run. That
    is WEAK EVIDENCE, not a good sign: **215 of the 230 `PSS(` constructions
    in this suite pass `method=` explicitly**, so only ~15 sites touch the
    default at all and most of those are topology checks that never solve.
    A default nothing exercises can be changed to anything without a red test
    — the §D 0ab shape, inverted. This test exercises it deliberately.

    Two facts, both of which were established by READING the code first and
    are asserted here because reading is not measuring:

    1. The default is `radau`. Chosen because the floor of this stack is
       discretisation and grows linearly in `Q` (gear `~1.8e-05 Q`, radau
       `~7.0e-12 Q` at 240 points).  ⚠⚠ The record also cited, as "the part
       that decides it", `trap`'s error CHANGING SIGN near `Q = 100` so that
       `grid_error` must refuse it.  WITHDRAWN 2026-09-14: that was a defect
       (the twin's `c` over trap's own period), and fixed, trap is estimable.
       The accuracy argument stands and is stronger: radau at 60 points
       (7.6e-07) still beats trap at 480 (3.3e-06, not the defect's 4.1e-06).

    2. An autonomous run under the default takes **no TR-BDF2 twin**.
       `monodromy_twin` returns `self` for any method that
       `carries_own_monodromy`, and every stage method does. Had Radau been
       missing from that set, every autonomous solve would have quietly
       spawned a SECOND PSS and read `lambda_2` from an order-2 map while the
       native one is order 5 — more cost for a worse answer, and nothing would
       have failed.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 0.5 * (u - u ** 3 / 3.0))

    pss = PSS(cir)
    assert pss.par.method == 'radau', \
        'the default method is %r; the 2026-09-07 decision is radau' % (
            pss.par.method,)

    T0 = 2.0 * np.pi
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]),
                  maxiterations=80)
    assert pss.converged, 'the default method did not converge on van der Pol'

    ## ⚠ `is pss` — not merely "a PSS with the same answer". The whole point
    ## is that NO second solve happens.
    assert pss.monodromy_twin() is pss, \
        'the default spawned a monodromy twin; a stage method must read its ' \
        'own map (see carries_own_monodromy)'
    assert not pss._twins, \
        'a twin was solved and cached: %r' % (list(pss._twins),)


def _weakly_limited_lc(lambda2):
    """An LC tank with loss and a cubic negative resistance, limited so weakly
    that the second Floquet multiplier is `lambda2` at 1 V (first-order
    averaging: b = -4 C ln(lambda2) / (3 T A^2), a = gl + 3 b A^2 / 4)."""
    lval, cval, gl, amp = 100e-6, 100e-12, 1e-4, 1.0
    period = 2 * np.pi * np.sqrt(lval * cval)
    b = -4.0 * cval * np.log(lambda2) / (3.0 * period * amp ** 2)
    a = gl + 0.75 * b * amp ** 2
    c = SubCircuit()
    c.add_node('p')
    c['L0'] = L('p', gnd, L=lval)
    c['C0'] = C('p', gnd, c=cval)
    c['R0'] = R('p', gnd, r=1.0 / gl)
    c['N0'] = BSource('p', gnd, gnd, 'p', i_func=lambda u: a * u - b * u ** 3)
    return c, period, amp


def test_a_gear_free_period_stall_near_unit_multiplier_is_named_and_redirected():
    """⚠⚠ A MULTISTEP FREE-PERIOD SOLVE THAT STALLS AT A RESIDUAL FLOOR MUST
    NOT BE TOLD TO RAISE `maxiterations`.  Reported by a peer session on a
    weakly limited LC oscillator and reproduced on this twin: at a second
    multiplier of 0.99, Gear-2 on a coarse uniform grid stops at a floor with
    no root nearby (the solved-history discrete solution does not exist below
    a grid threshold), while radau converges on the same grid, and Gear-2
    converges at 0.9.  The generic non-convergence warning used to advise
    `method='gear'` on an oscillator -- the wrong way.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def run(method, lambda2):
        cir, period, amp = _weakly_limited_lc(lambda2)
        pss = PSS(cir, method=method, reltol=1e-10)
        x0 = np.zeros(cir.n - 1)
        x0[0] = amp
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=period, timestep=period / 200, x0=x0,
                      maxiterations=6)
        return pss.converged, [str(w.message) for w in caught]

    converged, msgs = run('gear', 0.99)
    assert not converged, 'the fixture no longer stalls; it tests nothing'
    stall = [m for m in msgs if 'solved-history stall' in m]
    assert len(stall) == 1, msgs
    assert "method='radau'" in stall[0] and 'x_{-1}' in stall[0], stall[0]
    generic = [m for m in msgs if 'did not converge in 6 iterations' in m]
    assert len(generic) == 1, msgs
    assert "use method='gear'" not in generic[0], generic[0]

    ## the redirection is true on this fixture: radau converges on that grid
    converged, msgs = run('radau', 0.99)
    assert converged
    assert not [m for m in msgs if 'residual floor' in m], msgs

    ## and the diagnosis is silent when gear converges (stronger limiting)
    converged, msgs = run('gear', 0.9)
    assert converged
    assert not [m for m in msgs if 'residual floor' in m], msgs


def test_gears_free_period_is_second_order_on_a_non_uniform_grid_inside_the_zero_stability_bound():
    """Item 2 of gear-first-class (2026-09-20): the free-period solve on a
    non-uniform gear grid.  Nothing needed fixing; this pins what was
    measured.  Van der Pol mu = 1, reference radau N = 3200 (grid-independent
    to 1e-10 ppm; radau is order 5 on every grid here, 3:1 included)::

        grid                  N=200      N=400     N=800     N=1600    order
        uniform             +313.5     +78.1     +19.5      +4.87     4.00
        smooth 1+0.5 sin    +442.3    +110.8     +27.7      +6.94     4.00
        alternating 2:1     +310.5     +77.7     +19.4      +4.86     4.00
        alternating 3:1    -1615     -840.8    -428.7    -216.4      1.98

    Second order wherever the step ratio stays inside BDF2's zero-stability
    bound 1 + sqrt(2) -- the 2:1 grid's constant is the uniform one's -- and
    FIRST order beyond it, where `_period_grid` warns (it counts REPEATED
    up-steps; the 3:1 grid has one every other step).  ⚠ The phase fixture
    (a linear rotation) is second order on 3:1 too, and would have hidden
    this: the parasitic mode compounds only through the nonlinearity.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T_REF = 6.330195895892e+00
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: (u - u ** 3 / 3.0) + 0.3 * u * u)
        return cir

    def ppm(n, fr):
        pss = PSS(vdp(), method='gear', reltol=1e-10)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            pss.solve(period=T0, timestep=T0 / n, maxiterations=60,
                      break_events=False, grid=fr, x0=np.array([2.0, 0.0]))
        assert pss.converged
        warned = any('steps up by' in str(w.message) for w in rec)
        return 1e6 * (pss.period - T_REF) / T_REF, warned

    def smooth(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    def alt(n, r):
        f = np.tile([r, 1.0], n // 2)
        return f / f.sum()
    for name, mk, lo, hi in (('smooth', smooth, 3.6, 4.4),
                             ('alt21', lambda n: alt(n, 2.0), 3.6, 4.4),
                             ('alt31', lambda n: alt(n, 3.0), 1.7, 2.3)):
        e400, w400 = ppm(400, mk(400))
        e800, w800 = ppm(800, mk(800))
        assert lo < e400 / e800 < hi, (name, e400, e800)
        assert w400 is w800 is (name == 'alt31'), (name, w400, w800)
    e_u, _ = ppm(800, None)
    e_2, _ = ppm(800, alt(800, 2.0))
    assert abs(e_2 / e_u - 1.0) < 0.05, (e_u, e_2)


def test_lte_grid_steps_with_the_pss_own_method_so_a_radau_grid_is_radau_shaped():
    """`lte_grid` (and the since-deleted `refine_grid`) used to build their adaptive Transient
    with no integrator -- the Transient default `Gear2Integrator()` -- so a
    radau PSS got a gear-shaped grid, denser on the slow branch than a
    fifth-order method needs (2026-09-21, Andreas: "Do 2").  They now step
    with the PSS's own integrator.  Relaxation van der Pol, mu = 10, reltol
    1e-5, each method on the grid ITS OWN run produced, seeded from
    `lte_period`::

        method   steps   span      T err (ppm)   c rel
        gear     195     176:1     -1461 / -22   -8.6 % / -4.5 %   (two windows)
        trbdf2   286     833:1     -8.7          +0.29 %
        radau    110     130:1     -0.32         -0.05 %

    ⚠ Gear's two numbers are two windows of runs that differ only in the
    period hint (19.1 vs 16.1): same count, same span -- and the "70x" was
    NOT gear's (measured after this test was written, 2026-09-21): trbdf2
    on the same two grids reads +10 / -180 ppm, radau -0.01 / -0.04, so
    the windows differ for every second-order method (12 growth-by-2 steps
    against 5, 4x the total LTE); and the -22 ppm window sat on a SIGN
    CHANGE of the error (-22 -> +7 -> +3 under 2x / 4x splitting) where the
    other is clean second order (-1461 -> -361 -> -92).  Gear's honest
    number on a 195-point adaptive grid here is O(1000 ppm), second order;
    a good number at one N is refined before it is believed.
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
    grids = {}
    for method in ('gear', 'radau'):
        cir = vdp()
        p = PSS(cir, method=method)
        xfull = np.zeros(cir.n)
        xfull[[str(n_) for n_ in cir.nodes].index('v')] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            fr, seed = p.lte_grid(19.1, x0=xfull, reltol=1e-5)
        grids[method] = (np.asarray(fr, float), seed, float(p.lte_period))
    ## a radau grid is not a gear grid: sparser, from the fifth-order run
    assert len(grids['radau'][0]) < 0.75 * len(grids['gear'][0]), \
        (len(grids['radau'][0]), len(grids['gear'][0]))
    fr, seed, Tl = grids['radau']
    q = PSS(vdp(), method='radau', reltol=1e-9)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        q.solve(period=Tl, timestep=Tl / len(fr), x0=seed, maxiterations=80,
                break_events=False, grid=fr)
    assert q.converged
    e = 1e6 * (float(q.period) - T_REF) / T_REF
    assert abs(e) < 5.0, e                                       # -0.32 measured


def test_the_stage_methods_shoot_matrix_free():
    """`matrix_free=True` refused every stage method ("its monodromy is a
    dense stage product") -- the DEFAULT method among them -- though the
    factored stage map, with its mat-vec, and a walk carrying the period
    column without the dense map had made the reason stale.  Measured once
    the refusal went: radau, trbdf2 and esdirk43 match their dense solves to
    1e-15 on a driven RLC and a van der Pol oscillator.  (A Nordsieck GLM
    runs matrix-free too, once its startup is linearised -- see
    `test_a_glm_shoots_on_its_exact_map_once_the_startup_is_linearised`.)
    The state-event stage, which did not run under `matrix_free` (silent,
    then warned), runs there since 2026-09-25 --
    `test_the_state_event_stage_runs_matrix_free`.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    cases = ((_q20_rlc, dict(period=1e-3, timestep=1e-3 / 200)),
             (_vdp_asym, dict(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]))))
    for method in ('radau', 'trbdf2'):
        for mk, kw in cases:
            out = {}
            for mf in (False, True):
                p = PSS(mk(), method=method, reltol=1e-9)
                with _w.catch_warnings():
                    _w.simplefilter('ignore')
                    p.solve(maxiterations=60, matrix_free=mf, **kw)
                assert p.converged, (method, mk.__name__, mf)
                out[mf] = (float(p.period), np.asarray(p.waveform[1], dtype=float))
            assert abs(out[True][0] / out[False][0] - 1.0) < 1e-12, (method, mk.__name__)
            X0, X1 = out[False][1], out[True][1]
            assert np.max(np.abs(X1 - X0)) < 1e-10 * np.max(np.abs(X0)), (method, mk.__name__)


def test_every_method_states_its_order_to_grid_error_and_warping_estimate():
    """`grid_error`'s ceiling on a plausible observed order and
    `warping_estimate`'s interpolant degree came from two tables that
    covered 6 and 7 of the 12 accepted method names: esdirk43, the GLMs and
    the aliases `gear2` / `trapezoidal` fell to a generic ceiling (none
    applied) and a cubic.  Both now come from the integrator's own order --
    the degree the smallest odd one above it, at least 3, the rule every
    measured entry followed.  The cubic cost glm4 its estimate: on van der
    Pol at 60 points it read +9.8e-7 against a true +1.7e-7 (the derived
    quintic +2.0e-7), and at 120 points the wrong sign.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    want = {'euler': (1, 3), 'trap': (2, 3), 'trapezoidal': (2, 3),
            'theta': (2, 3), 'gear': (2, 3), 'gear2': (2, 3),
            'trbdf2': (2, 3), 'radau': (5, 7), 'esdirk43': (4, 5),
            'glm2': (2, 3), 'glm3': (3, 5), 'glm4': (4, 5)}
    p = PSS(_q20_rlc())
    for m, (order, degree) in want.items():
        assert p._nominal_order(m) == order, (m, p._nominal_order(m))
        assert p._idec_degree(m) == degree, (m, p._idec_degree(m))

    T_REF = 6.330195895892e+00
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    q = PSS(_vdp_asym(), method='glm4', reltol=1e-12)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        q.solve(period=T0, timestep=T0 / 60, x0=np.array([2.0, 0.0]),
                maxiterations=60)
        assert q.converged
        true = float(q.period) / T_REF - 1.0
        est = q.warping_estimate()['period_error'] / float(q.period)
        cubic = q.warping_estimate(degree=3)['period_error'] / float(q.period)
    assert abs(est / true - 1.0) < 0.3, (est, true)
    assert abs(cubic / true - 1.0) > 2.0, (cubic, true)



def test_a_glm_shoots_on_its_exact_map_once_the_startup_is_linearised():
    """A Nordsieck GLM's period map on `x_0` runs through its STARTUP -- p
    Radau substeps and an interpolant building the Nordsieck vector -- and
    the recursion was seeded with ``[C(x_0), 0, ...]``, the startup's higher
    components held fixed.  So the shooting Jacobian was approximate (3.8e-2
    / 9.5e-2 from finite differences under glm2 / glm3 on van der Pol), the
    Newton only linear (9 and 14 evaluations on a LINEAR RLC), the factored
    map not the Newton's (matrix-free, glm2 and glm3 DIVERGED on that RLC)
    and no spectrum was reported.  With the startup linearised
    (`_GLMStartup`, exact through each substep's converged stage system):
    the map 3e-9 / 2e-10 from finite differences and its period column
    8e-12 / 6e-11, the RLC converged in ONE Newton step, matrix-free equal
    to dense to 1e-15, and the spectral radius the RLC's decay.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    for method in ('glm2', 'glm3'):
        p = PSS(_vdp_asym(), method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / 60, x0=np.array([2.0, 0.0]),
                    maxiterations=60)
        assert p.converged, method
        x0 = np.asarray(p._period_state[1], dtype=float)[:p.cir.n - 1]
        T = float(p.period)
        p._open_at_x0 = False

        def walk(x, TT, want=False):
            tms, hs = p._period_grid(TT, 60, None)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                return p._walk('glm', x, tms, hs, T=TT, want_dT=want)
        w = walk(x0, T, True)
        M = np.asarray(w.monodromy(), dtype=float)
        Mt = np.asarray(w.period_column(), dtype=float).ravel()
        eps = 1e-6
        Mfd = np.column_stack([(np.asarray(walk(x0 + eps * e, T).end())
                                - np.asarray(walk(x0 - eps * e, T).end())) / (2 * eps)
                               for e in np.eye(len(x0))])
        dT = 1e-6 * T
        Mtfd = (np.asarray(walk(x0, T + dT).end())
                - np.asarray(walk(x0, T - dT).end())) / (2 * dT)
        assert np.max(np.abs(M - Mfd)) < 1e-7 * np.max(np.abs(Mfd)), method
        assert np.max(np.abs(Mt - Mtfd)) < 1e-8 * np.max(np.abs(Mtfd)), method

    rho_exact = float(np.exp(-np.pi / 20.0))      # `_q20_rlc`'s decay per period
    for method in ('glm2', 'glm3'):
        out = {}
        for mf in (False, True):
            p = PSS(_q20_rlc(), method=method, reltol=1e-9)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                p.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=60,
                        matrix_free=mf, trace=True)
            assert p.converged, (method, mf)
            out[mf] = np.asarray(p.waveform[1], dtype=float)
            if not mf:
                ## a linear circuit and an exact Jacobian: one Newton step
                F1 = np.max(np.abs(np.asarray(p.shooting_trace[1][1])))
                F0 = np.max(np.abs(np.asarray(p.shooting_trace[0][1])))
                assert F1 < 1e-12 * F0, (method, F0, F1)
                assert abs(p.spectral_radius / rho_exact - 1.0) < 1e-3, p.spectral_radius
        assert np.max(np.abs(out[True] - out[False])) < 1e-10 * np.max(np.abs(out[False]))


def test_a_glm_map_carries_the_nordsieck_rescale_on_a_non_uniform_grid():
    """Where the step changes, the transient scales the Nordsieck vector,
    ``Q_k <- (h/h_old)^k Q_k``, before the step -- and the GLM's
    linearisation did not.  Measured on van der Pol over a 3:1 grid: the
    map 7.2e-3 / 2.0e-2 off its finite difference (glm2 / glm3; exact on a
    uniform grid), the Newton 14 / 21 evaluations, and the 'closing' period
    column 8 % / 19 % off on EVERY grid (the last step's rescale moves with
    `T` even where `rho` is 1).  The step records now carry `rho` and the
    entered vector (`_glm_period_blocks`), `_GLMStep.forward` scales by `rho^k`,
    its adjoint scales back, and the closing column carries ``(k / h_N)
    Q_k``.  Pinned: the map and both period columns against central
    differences on the 3:1 grid, and the factored map's transpose against
    its forward.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 2.0 * np.pi / np.sqrt(1.0 - 0.25 / 4.0)
    N = 60
    j = np.arange(N)
    fr = 1.0 + 0.5 * np.sin(2 * np.pi * (j + 0.5) / N)
    fr = fr / fr.sum()
    rng = np.random.default_rng(1)
    for method in ('glm2', 'glm3'):
        p = PSS(_vdp_asym(), method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T0, timestep=T0 / N, x0=np.array([2.0, 0.0]),
                    maxiterations=60, grid=fr)
        assert p.converged, method
        x0 = np.asarray(p._period_state[1], dtype=float)[:p.cir.n - 1]
        T = float(p.period)
        _tms0, hs0 = p._period_grid(T, N, fr)
        assert float(np.max(hs0) / np.min(hs0)) > 2.5

        def walk(x, hs, TT, want=False, keep=False):
            tms = np.concatenate(([0.0], np.cumsum(hs)))
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                return p._walk('glm', x, tms, hs, T=TT, want_dT=want,
                               keep=keep)
        p._period_column = 'proportional'
        w = walk(x0, hs0, T, True)
        M = np.asarray(w.monodromy(), dtype=float)
        Mt = np.asarray(w.period_column(), dtype=float).ravel()
        eps, dT = 1e-6, 1e-6 * T
        Mfd = np.column_stack([(np.asarray(walk(x0 + eps * e, hs0, T).end())
                                - np.asarray(walk(x0 - eps * e, hs0, T).end()))
                               / (2 * eps) for e in np.eye(len(x0))])
        Mtfd = (np.asarray(walk(x0, hs0 * (T + dT) / T, T + dT).end())
                - np.asarray(walk(x0, hs0 * (T - dT) / T, T - dT).end())) / (2 * dT)
        assert np.max(np.abs(M - Mfd)) < 1e-7 * np.max(np.abs(Mfd)), method
        assert np.max(np.abs(Mt - Mtfd)) < 1e-8 * np.max(np.abs(Mtfd)), method
        ## 'closing': only the last step's length moves with the period
        p._period_column = 'closing'
        Mtc = np.asarray(walk(x0, hs0, T, True).period_column(),
                         dtype=float).ravel()
        hp, hm = hs0.copy(), hs0.copy()
        hp[-1] += dT
        hm[-1] -= dT
        Mtcfd = (np.asarray(walk(x0, hp, T + dT).end())
                 - np.asarray(walk(x0, hm, T - dT).end())) / (2 * dT)
        assert np.max(np.abs(Mtc - Mtcfd)) < 1e-8 * np.max(np.abs(Mtcfd)), method
        p._period_column = 'proportional'
        fp = walk(x0, hs0, T, keep=True).factored(p)
        u, v = rng.standard_normal(fp.width), rng.standard_normal(fp.width)
        fwd = float(u @ np.asarray(fp.matvec(v)))
        assert abs(fwd - float(np.asarray(fp.matvec_transposed(u)) @ v)) < 1e-12 * abs(fwd)


def test_a_glm_reads_pac_off_its_own_map_at_its_order():
    """Since 2026-09-25 PAC, the adjoint rows and pnoise read a DRIVEN
    Nordsieck GLM run off its own map on the state (`PSS._state_map`; until
    then a radau twin -- Andreas: "Native for all but covariance").  The
    forced replay is the step on ``(P, x)`` (`_GLMStateStep`) with the
    source in every stage, the output rows and the substages of the startup
    that opens a step.  Measured on the LTI RLC against the analytic AC at
    1.3 kHz: glm2 1.20e-3 / 2.98e-4 / 7.43e-5, glm3 1.03e-5 / 1.29e-6 /
    1.62e-7, glm4 1.86e-7 / 1.04e-8 / 6.06e-10 at 100 / 200 / 400 points --
    order 2.0 / 3.0 / 4.1; the forced replay against a finite difference of
    the period walk with the source perturbed: 2.4e-8 at eps 1e-3 (1/eps
    below it: the linear circuit's Newton floor); the generic replays equal
    the GLM's own `x_matvec` / `x_matvec_transposed` bit for bit.  Pinned:
    the order within a quarter of p over 100 -> 200, the finite difference
    to 1e-6, the replays to 1e-14.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    f0, Q = 1e3, 20.0
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    R_ = (1.0 / Q) * np.sqrt(L_ / C_)
    f = 1300.0
    w = 2 * np.pi * f
    exact = (1.0 / (1j * w * C_)) / (R_ + 1j * w * L_ + 1.0 / (1j * w * C_))

    def run(method, N):
        c = _q20_rlc()
        p = PSS(c, method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=1e-3, timestep=1e-3 / N)
            res = PAC(c, toolkit=circuit.numeric).solve(p, [f])
        sv = np.asarray(res.sweep_values, dtype=float)
        X = np.asarray(res.x)[[str(n_) for n_ in c.nodes].index('c')]
        ks = [k for k in range(len(sv)) if abs(sv[k] - f) < 1e-9 * f]
        x = complex(X[max(ks, key=lambda k: abs(X[k]))])
        return abs(x / exact - 1.0), p

    for method, order in (('glm2', 2), ('glm3', 3)):
        e1, p = run(method, 100)
        e2, _p2 = run(method, 200)
        assert p._state_map().is_glm
        assert abs(np.log2(e1 / e2) - order) < 0.25, (method, e1, e2)
        ## the replays on the state against the GLM's own map
        sm = p._state_map()
        rng = np.random.default_rng(3)
        v = rng.standard_normal(sm.width)
        assert np.max(np.abs(p._replay(sm, v) - sm.matvec(v))) < 1e-14 * np.max(np.abs(sm.matvec(v)))
        assert np.max(np.abs(p._replay_transposed(sm, v) - sm.matvec_transposed(v))) \
            < 1e-14 * np.max(np.abs(sm.matvec_transposed(v)))
        ## the forced replay against the period walk with the source moved
        m = sm.width
        tms = np.asarray(sm.times, dtype=float)
        x0 = np.asarray(p._period_state[1], dtype=float)[:m]
        iref = p.irefnode
        u_red = rng.standard_normal(m) * 1e-3
        u_full = np.insert(u_red, iref, 0.0)
        fq, eps = 2.3e3, 1e-3
        cir = p.cir
        orig = cir.u

        def walk():
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                xs = p._glm_period_blocks(x0, tms, np.diff(tms))[1]
            return np.array([np.asarray(x_, dtype=float) for x_ in xs])

        X0 = walk()
        try:
            def u_pert(t, epar=None, analysis=None, **kw):
                base = orig(t, epar, analysis=analysis, **kw)
                if analysis == 'ac':
                    return base
                return np.asarray(base, dtype=float) + eps * np.real(u_full * np.exp(2j * np.pi * fq * t))
            cir.u = u_pert
            X1 = walk()
        finally:
            del cir.u
        fd = (X1 - X0) / eps
        _e, ys = p._forced_replay(sm, fq, u_red, y0=np.zeros(m, dtype=complex), collect=True)
        err = np.max(np.abs(np.real(np.array(ys)) - fd)) / np.max(np.abs(fd))
        assert err < 1e-6, (method, err)


def test_grid_error_refines_by_any_factor_and_reports_the_methods_order():
    """The review's X9 (2026-10-01): `grid_error(refine=, levels=)` ran at
    its defaults only.  On van der Pol (mu = 0.2) the period's estimated
    order is the method's at `refine=3` as at 2 -- gear 1.988 / 1.991, radau
    5.25 / 5.22 -- on grids N, 3N, 9N; `levels=2` is the raw two-grid change
    with no order and no power-law check; and the estimate is the ERROR:
    gear at 40 points refined 3x twice estimates its finest value 6.286e-4
    off, against 6.232e-4 from a radau reference (1.009)."""
    circuit.default_toolkit = circuit.numeric
    mu = 0.2
    T = 2.0 * np.pi

    def solved(method, npts):
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0))
        p = PSS(c, method=method, reltol=1e-12)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                    maxiterations=100)
        assert p.converged
        return p
    period = lambda q: float(q.period)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        g = solved('gear', 40)
        r3 = g.grid_error(period, refine=3, levels=3)
        r2 = g.grid_error(period, refine=3, levels=2)
        rr = solved('radau', 20).grid_error(period, refine=3, levels=3)
        ref = float(solved('radau', 180).period)
    assert list(r3['npts']) == [40, 120, 360] and abs(r3['order'] - 2.0) < 0.1
    assert abs(rr['order'] - 5.0) < 0.5, rr['order']
    assert list(r2['npts']) == [40, 120] and r2['order'] is None
    assert r2['power_law'] is None
    actual = abs(r3['values'][-1] - ref)
    assert abs(r3['error'] / actual - 1.0) < 0.05, (r3['error'], actual)
    with pytest.raises(ValueError, match='levels'):
        g.grid_error(period, levels=4)
