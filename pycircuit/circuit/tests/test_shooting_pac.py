"""Shooting tests: shooting pac.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    ConvergenceWarning,
    UsageWarning,
)
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
import warnings
from pycircuit.circuit.tests._shooting_elements import _SwitchHdl
from pycircuit.circuit.tests._shooting_fixtures import (_adjoint_ladder,
    _comparator_relaxation_oscillator,
    _cos2_mixer,
    _diode_mixer,
    _exact_relaxation_oscillator_model,
    _exact_relaxation_oscillator_ppv,
    _lc_osc,
    _pac_circuit,
    _pwm_loop,
    _rc_ladder,
    _relaxation_oscillator_seed,
    _solve_slow)


def test_PAC_runs_on_a_nonlinear_circuit_without_forming_the_operator():
    """⚠ PAC IS NO LONGER WITHDRAWN — this test used to assert that it was.

    Stage 11 withdrew it for forming the whole `(N m) x (N m)` operator
    densely: 419.5 GiB at `N = 137`, `m = 1000`. The withdrawal note said
    the body was "the starting point for a matrix-free rewrite", and that
    turned out to be exactly right — the operator in it was CORRECT, just
    un-preconditioned. See `PAC`'s docstring and
    `test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`.

    This keeps the withdrawal's own circuit — a diode, so genuinely
    nonlinear and genuinely time-varying, which is the case PAC exists for
    and the one the AC gate cannot cover — and asserts the two things the
    rewrite promises: it runs, and it never allocates the operator.

    ⚠ THE MEMORY ASSERTION IS THE POINT, not the numbers. The dense route
    would allocate `(N m)^2` complex entries; this circuit is small enough
    that doing so would succeed and pass a value check silently. So the
    test counts what is stored: `N` factorisations of an `m x m` matrix,
    which is `O(N m^2)` and not `O(N^2 m^2)`.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    pss = PSS(cir, method='gear')
    with quiet(AccuracyWarning):
        pss.solve(period=1e-6, timestep=1e-6 / 40)
    assert pss.converged

    pac = PAC(cir)
    with quiet():
        res = pac.solve(pss, np.array([1e6, 3e6]))
    X = np.asarray(res.x, dtype=complex)
    assert np.isfinite(X).all(), 'PAC returned non-finite entries'
    assert np.abs(X).max() > 0, \
        'PAC returned all zeros on a driven circuit -- the source vector is ' \
        'probably being read positionally again (see the AC gate)'

    ## the operator is never formed: what is stored is per-STEP, m x m
    fp = pss.factored_period()
    m = cir.n - 1
    N = len(fp.steps)
    assert N > 1 and fp.width <= 2 * m, \
        'the factored period should hold N per-step m x m factorisations, ' \
        'not one (N m) x (N m) matrix'
    assert res.info['matvecs'] < 4 * fp.width, \
        'PAC used %d matvecs for 2 frequencies on an m=%d circuit; forming ' \
        'the monodromy would take %d, and the point is not to' \
        % (res.info['matvecs'], m, fp.width)


@unittest.skip("Superseded by the PAC gate tests at the end of this file")
def test_PAC():
    circuit.default_toolkit = circuit.numeric
    N = 10
    fc = 1e6

    cir = SubCircuit()
    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=fc, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2,gnd)
    cir['C'] = C(2,gnd, c=1e-12)
    
    pss = PSS(cir)

    res = pss.solve(period=1/fc, timestep = 1/(fc*N))

    pac = PAC(cir)
    res = pac.solve(pss, freqs = fc + np.array([1e3, 2e3, 4e3]))
    
    assert False, "Test should compare with a reference simulation"


def test_a_failed_inner_krylov_solve_says_so():
    """A Krylov failure must name itself, not become 'No convergence'.

    ⚠ THE INNER SOLVE'S VERDICT USED TO BE DISCARDED.  `gmres` returns an
    `info` saying whether it converged, broke down, or ran out of cycles,
    and dropping it meant a matrix-free run could spend its whole budget --
    scipy's `maxiter` counts RESTART CYCLES, so the pair multiplies, and the
    earlier `restart=200, maxiter=400` was a worst case near 80 000 matvecs,
    each a full replay of the period -- and then report the same generic
    outer message as a circuit that simply needed another iteration.

    This file's standard is that a failure says what happened, which is why
    `T = 0` and the singular free-period Jacobian are named.  The inner
    solve is now held to it too.  Starved of cycles on purpose, because the
    honest way to test a failure path is to cause the failure.
    """
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_rc_ladder(20), method='gear', reltol=1e-6)
    old_cycles, old_restart = pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART
    try:
        ## one cycle of a one-dimensional Krylov space cannot solve this
        pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART = 1, 1
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=1e-3, timestep=1e-3 / 25, maxiterations=5,
                      matrix_free=True)
    finally:
        pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART = old_cycles, old_restart

    assert not pss.converged, \
        'a starved inner solve still reported convergence'
    msgs = [str(w.message) for w in caught]
    assert any('inner solve did not converge' in m for m in msgs), \
        'the Krylov failure was not named; warnings were %r' % msgs
    assert any('matrix_free=False' in m for m in msgs), \
        'the message does not say what to do instead'


def _pac_y0(pss, freq, per):
    """`y_0` by a DENSE solve of the same `m x m` system PAC solves.

    Dense on purpose: this is the reference the matrix-free path is
    measured against, so it must not share its solver.
    """
    fp = pss.factored_period()
    irn = pss.irefnode
    (u_ac,) = remove_row_col((pss.cir.u(0, analysis='ac'),), irn,
                             pss.toolkit)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    w, _ = pss._forced_replay(fp, freq, u_ac)
    alpha = np.exp(-2j * np.pi * freq * per)
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    return np.linalg.solve(np.eye(n) - alpha * M, alpha * np.asarray(w)), M


def _pac_vs_ac(method, npts, freq=700.0, per=1e-3, x0_unknown=False):
    from pycircuit.circuit.analysis_ss import AC
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / npts, maxiterations=40,
                  x0_unknown=x0_unknown)
    assert pss.converged, '%s at %d points did not converge' % (method, npts)
    y0, _M = _pac_y0(pss, freq, per)
    irn = pss.irefnode
    xac = np.asarray(AC(cir, toolkit=circuit.numeric).solve(freqs=freq).x,
                     dtype=complex).ravel()
    xac = np.concatenate((xac[:irn], xac[irn + 1:]))
    return (np.linalg.norm(np.asarray(y0)[:m] - xac) / np.linalg.norm(xac),
            pss)


def test_pac_agrees_with_the_ac_analysis_on_a_linear_circuit():
    """⚠ THE GATE FOR PAC, against a reference it cannot influence.

    On a LINEAR circuit the periodic operating point is irrelevant to the
    small signal: the linearisation is constant, every sideband but `l = 0`
    vanishes, and the LPTV response collapses to the LTI transfer function.
    So `AC` — a different analysis, a different code path, an `(sC + G)`
    solve that never touches a monodromy — is the answer PAC must produce.

    This is the check the withdrawn implementation never had: its only test
    was `@unittest.skip('Skip failing test')`.

    ⚠ AND IT IS WHAT CAUGHT THE DEAD BODY'S REAL DEFECT. That body read its
    source vector as `self.cir.u(0, analysis_name)` — POSITIONALLY, into a
    signature whose second parameter is `epar`. It took the transient source
    at `t = 0`, which is zero for every sinusoid, so the whole analysis
    would have returned zeros. Reproduced here before it was fixed: `|PAC|`
    exactly 0 against `|AC|` of 1.01, at every frequency and every method.
    """
    rel, _pss = _pac_vs_ac('gear', 1000)
    assert rel < 1e-4, \
        'PAC disagrees with the AC analysis by %.3e on a LINEAR circuit, ' \
        'where the two are solving the same problem by different routes' % rel
    ## ⚠ PRECONDITION / BLIND TO (2026-09-05): AC is the right reference
    ## only because this circuit does NOT convert -- the l=1 sideband must
    ## be ~0 here, and asserting it says so.  It also means this gate
    ## CANNOT see the sideband decomposition; the converting-circuit gates
    ## that do are `test_pac_reports_sidebands_at_the_right_frequencies_`
    ## `and_conjugates_the_fold` and `test_a_switched_capacitor_holds_kTC_`
    ## `with_per_step_CY`.
    _pac = PAC(_pss.cir, toolkit=circuit.numeric)
    _oc = [str(n) for n in _pss.cir.nodes].index('c')
    _H = np.asarray(_pac.adjoint_sideband_row(_pss, 700.0, _oc,
                                              sidebands=[0, 1])[0])
    assert np.linalg.norm(_H[1]) < 1e-6 * np.linalg.norm(_H[0]), \
        'the l=1 sideband is %.3e of l=0 on a circuit that must not ' \
        'convert; AC is not the right reference if it does' \
        % (np.linalg.norm(_H[1]) / np.linalg.norm(_H[0]))


@pytest.mark.parametrize('method,tol', [('trbdf2', 1e-4), ('radau', 1e-9)])
def test_pac_forward_replay_works_over_the_stage_methods(method, tol):
    """The FORWARD driven replay (PAC.solve) runs over the self-starting stage
    methods and agrees with the AC analysis on a linear circuit.

    The forward path was the last TR-BDF2/Radau deferral: a stage step injects
    the source at THREE abscissae (``A (x) B``), which the LMM one-injection
    fold cannot carry.  ``_forced_replay_{trbdf2,radau}`` carry it -- the exact
    transpose of the (already-shipped) adjoint fold.  On a linear circuit PAC
    must reduce to AC, a reference the shooting path cannot influence; Radau
    (order 5) reaches it far tighter than TR-BDF2 (order 2) at the same grid.
    """
    with quiet():
        rel, _pss = _pac_vs_ac(method, 200)
    assert rel < tol, \
        '%s PAC disagrees with AC by %.3e on a LINEAR circuit' % (method, rel)


@pytest.mark.parametrize('method,x0_unknown,expect', [
    ('trap', False, 'first'),
    ('trap', True, 'second'),
    ('euler', False, 'first'),
    ('euler', True, 'first'),
])
def test_pac_order_is_lost_to_the_manufacturing_step(method, x0_unknown,
                                                     expect):
    """⚠ PAC ON THE PLAIN PATH IS FIRST ORDER WHATEVER THE METHOD.

    The plain kept walk takes one step OUTSIDE its loop to
    manufacture a history, and that step is not in `steps`. For the
    homogeneous map its effect is folded into the `opening` triple — the
    documented flat-history approximation. For the DRIVEN map it also means
    THE SOURCE IS NEVER APPLIED THERE: one step of `u` out of `N`, a
    relative O(h), which drags a second-order method down to first.

    ⚠ THE RATE IS THE ASSERTION, NOT THE SIZE. A constant-factor error is
    invisible to "is it small" and unmissable to "does it fall" — the same
    reason the saltation falsifier and the A&T period column are rate
    checks. Measured against the AC analysis at 700 Hz, per doubling:

        trap, plain            2.00x   4.13e-03 at 250 points
        trap, x0_unknown=True  4.00x   1.09e-04 at 250 points
        euler, either          2.00x   1.40e-02, identical to five digits

    ⚠ EULER IS THE CONTROL AND IT IS WHY THIS IS THE MANUFACTURING STEP.
    `x0_unknown` does not move euler's rate at all — it is a first-order
    method either way — so the trapezoidal gain cannot be some other thing
    the formulation does. The trajectory is not implicated either: trap's
    WAVEFORM converges at ~4.2x per doubling on this circuit with or
    without the manufacturing step.
    """
    ## ⚠ BLIND TO (2026-09-05): the linear resonator of `_pac_circuit`
    ## does not convert, so this gate measures ORDER, never the sideband
    ## decomposition -- for that see the converting-circuit gates named in
    ## `test_pac_agrees_with_the_ac_analysis_on_a_linear_circuit`.
    r1, _ = _pac_vs_ac(method, 250, x0_unknown=x0_unknown)
    r2, _ = _pac_vs_ac(method, 500, x0_unknown=x0_unknown)
    ratio = r1 / r2
    if expect == 'second':
        assert 3.4 < ratio < 4.6, \
            '%s/x0_unknown=%s: %.4e -> %.4e is %.2fx per doubling, not the ' \
            "method's own second order. If the manufacturing step now " \
            'carries its source this test should be REWRITTEN, not relaxed' \
            % (method, x0_unknown, r1, r2, ratio)
    else:
        assert 1.7 < ratio < 2.4, \
            '%s/x0_unknown=%s: %.4e -> %.4e is %.2fx per doubling, not the ' \
            'first order the dropped manufacturing source predicts' \
            % (method, x0_unknown, r1, r2, ratio)


def test_pac_recycling_matches_the_per_frequency_solve():
    """ONE Krylov subspace for the whole sweep — Telichevesky et al. Thm 1.

    `A(alpha) = I - alpha M`, so `span{r, Ar, A^2 r, …}` is the Krylov
    space of `M` and does not depend on `alpha` at all. A basis built once
    therefore serves every frequency in the sweep, and each frequency costs
    a small dense least-squares over it rather than its own run of
    full-period replays.

    ⚠ THE RIGHT-HAND SIDE IS NOT SHARED, and that is why the implementation
    checks rather than assumes: `w(f)` genuinely differs per frequency, so
    the shared span is not the space GMRES would have picked for any but
    the first. It minimises the TRUE residual over the span and extends the
    basis — one matvec, kept for every later frequency — until every
    frequency is inside tolerance. The answer is therefore never worse than
    the per-frequency solve; only the count moves. Measured on RC ladders,
    24 frequencies:

        m=4    72 matvecs -> 6    12.0x
        m=14  168 matvecs -> 15   11.2x
        m=42  302 matvecs -> 26   11.6x

    with the two routes agreeing to 1.2e-13.
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, vac=1.0, freq=1e3)
    for k in range(6):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-7)
    per = 1e-3
    pss = PSS(c, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 150, maxiterations=40)
    assert pss.converged

    freqs = np.logspace(1, 4, 12)
    out = {}
    for rec in (False, True):
        pac = PAC(c, toolkit=circuit.numeric)
        with quiet():
            res = pac.solve(pss, freqs, recycle=rec)
        out[rec] = (np.asarray(res.x, dtype=complex), res.info['matvecs'])

    a, na = out[False]
    b, nb = out[True]
    rel = np.linalg.norm(a - b) / np.linalg.norm(a)
    assert rel < 1e-9, \
        'recycled sweep disagrees with the per-frequency solve by %.3e -- ' \
        'the shared subspace is supposed to change the COST and not the ' \
        'answer' % rel
    assert nb < na, \
        'recycling used %d matvecs against %d for the per-frequency solve, ' \
        'so it is not recycling anything' % (nb, na)


def test_pac_refuses_what_it_cannot_answer():
    """The two ways to ask PAC a question that has no answer.

    Both were reachable in the withdrawn body, and neither announced
    itself: a source with no `vac` gave a zero right-hand side and returned
    zeros, and nothing checked that the operating point had converged.
    """
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    ## (a) no operating point yet
    cir = _pac_circuit()
    pss = PSS(cir, method='gear', reltol=1e-10)
    with pytest.raises(RuntimeError, match='call solve'):
        pss.factored_period()

    ## (b) an operating point, but no small-signal source
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * 1e3) ** 2 * 1e-3)
    silent = SubCircuit()
    silent.add_node('a'); silent.add_node('b')
    silent['vs'] = VSin('a', gnd, va=1.0, vac=0.0, freq=1e3)
    silent['R'] = R('a', 'b', r=(1.0 / 20.0) * np.sqrt(L_ / C_))
    silent['L'] = L('b', 'c', L=L_)
    silent['C'] = C('c', gnd, c=C_)
    pss2 = PSS(silent, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss2.solve(period=per, timestep=per / 150, maxiterations=40)
    assert pss2.converged
    with pytest.raises(ValueError, match='identically zero'):
        PAC(silent, toolkit=circuit.numeric).solve(pss2, [700.0])


def test_the_adjoint_row_is_m_forward_solves_in_one():
    """⚠ THE MANY-TO-ONE IDENTITY, which is what makes pnoise affordable.

    pnoise is hundreds of sources into one output. Forward, that is one
    solve PER SOURCE, and recycling does not help because the right-hand
    side is what changes. Okumura et al. (1993) reach for the adjoint for
    exactly this reason -- "it is efficient to use the adjoint method …
    because circuits have many noise sources" -- and

        d^T y_0 = alpha ((I - alpha M)^-T d)^T W u

    turns the whole row into ONE transposed solve plus a reverse replay.

    ⚠ THE ASSERTION IS THE COUNT, NOT THE CLOCK. Measured speedup grows
    linearly with `m` (7.0x / 16.9x / 40.0x at m = 6 / 14 / 32, agreement
    ~1e-15), which is the right shape -- but this machine runs more than
    one agent and a wall-clock ratio here would be measuring the neighbours.
    So the test compares matvec counts, which no concurrent load can move.

    ⚠ AND IT COSTS NO NEW MACHINERY. The kept walk (`keep=True`) stores
    every step's factorisation and every factorisation already solves
    transposed, so the reverse pass needs no reverse integrator -- the
    thing Demir & Roychowdhury call "often unavailable even in existing
    time-domain simulators".
    """
    import scipy.sparse.linalg as spla
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder()
    per, f = 1e-3, 700.0
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 120, maxiterations=40)
    assert pss.converged
    fp = pss.factored_period()
    n = fp.width
    alpha = np.exp(-2j * np.pi * f * per)

    pac = PAC(cir, toolkit=circuit.numeric)
    row, row_info = pac.adjoint_transfer_row(pss, f, 1)
    adjoint_matvecs = row_info['matvecs']

    ## the forward route: one solve per source
    d = np.zeros(n)
    d[1] = 1.0
    fwd = np.zeros(m, dtype=complex)
    fwd_matvecs = [0]

    def _mv(v):
        fwd_matvecs[0] += 1
        return np.asarray(v) - alpha * fp.matvec(v)

    for i in range(m):
        u = np.zeros(m)
        u[i] = 1.0
        w, _ = pss._forced_replay(fp, f, u)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        y0, info = spla.gmres(A, alpha * np.asarray(w), rtol=1e-13,
                              restart=min(n, 50), maxiter=min(n, 200))
        assert info == 0
        fwd[i] = d @ y0

    rel = np.linalg.norm(row - fwd) / np.linalg.norm(fwd)
    assert rel < 1e-11, \
        'the adjoint row disagrees with %d forward solves by %.3e -- one ' \
        'of the two transposes is wrong, and the forward route is the one ' \
        'with an independent reference behind it' % (m, rel)
    assert adjoint_matvecs < fwd_matvecs[0], \
        'the adjoint used %d matvecs against %d for the forward route at ' \
        'm=%d; the whole point is that it does not scale with the number ' \
        'of sources' % (adjoint_matvecs, fwd_matvecs[0], m)


def test_the_adjoint_row_no_longer_refuses_the_plain_path():
    """⚠ THIS TEST USED TO ASSERT THE OPPOSITE, and the inversion is the
    record rather than a slip.

    It read: *"It needs the reverse replay, which is Gear-2 only"*, and
    demanded a `NotImplementedError` matching `solved-history`. That was
    true when written. B8 gave the one-step companions their own reverse
    recursion and then WIRED IT THROUGH, so the refusal is gone and the
    row must now compute -- on the same fixture, under the same method,
    where it previously raised.

    Kept as an inverted assertion rather than deleted: a deleted test
    leaves no evidence the capability was ever absent, and this one
    dates the change.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    per = 1e-3
    pss = PSS(cir, method='trap', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 120, maxiterations=40)
    assert pss.converged
    row = np.asarray(PAC(cir, toolkit=circuit.numeric)
                     .adjoint_transfer_row(pss, 700.0, 1)[0])
    assert row.shape == (cir.n - 1,), \
        'the adjoint row has the wrong width on the plain path: %r' \
        % (row.shape,)
    assert np.all(np.isfinite(row)), 'the plain adjoint row is not finite'
    assert float(np.max(np.abs(row))) > 0.0, \
        'the plain adjoint row came back all zeros, which is what a ' \
        'silently skipped reverse pass would produce'


def _sideband_forward(pss, fp, freq, k, l, N, alpha, A):
    """`H_l` for every source, the expensive way: one driven solve each."""
    m = pss.cir.n - 1
    tms = np.asarray(fp.times, dtype=float)
    w0 = 2.0 * np.pi / float(fp.T)
    out = np.zeros(m, dtype=complex)
    for i in range(m):
        u = np.zeros(m)
        u[i] = 1.0
        w, _ = pss._forced_replay(fp, freq, u)
        y0 = np.linalg.solve(A, alpha * np.asarray(w))
        _e, ys = pss._forced_replay(fp, freq, u, y0=y0, collect=True)
        y = [np.asarray(y0)[:m]] + [np.asarray(v)[:m] for v in ys]
        ## the DFT of `v = y exp(-j w t)`, which is what is T-periodic
        out[i] = sum(np.exp(-1j * (l * w0 + 2.0 * np.pi * freq) * tms[j])
                     * y[j][k] for j in range(N)) / N
    return out


def _sideband_pair(cir, per, freq, k, ls, npts=60, reltol=1e-11):
    """The adjoint rows and the `m` forward driven solves they replace."""
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method='gear', reltol=reltol)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    fp = pss.factored_period()
    n = fp.width
    N = len(fp.steps)
    alpha = np.exp(-2j * np.pi * freq * per)
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    A = np.eye(n) - alpha * M
    rows = PAC(cir, toolkit=circuit.numeric).adjoint_sideband_row(
        pss, freq, k, ls)[0]
    fwds = np.array([_sideband_forward(pss, fp, freq, k, l, N, alpha, A)
                     for l in np.atleast_1d(ls)])
    return pss, rows, fwds


@pytest.mark.parametrize('l', [0, 1, -2])
def test_the_sideband_adjoint_needs_both_terms(l):
    """⚠ A SIDEBAND IS A FUNCTIONAL OVER THE WHOLE PERIOD, not a value at
    one instant, and that changes the adjoint in two ways.

    `H_l` is the `l`-th Fourier coefficient of `v(t) = y(t) exp(-j w t)`,
    the part that is actually T-periodic. Every state on the trajectory
    contributes, so the reverse pass takes an injection at EVERY step
    rather than a seed at the end — and the initial state `y_0` is itself a
    function of the source through the periodic boundary condition, which
    is a SECOND term:

        dH/du = [forced part, from the injected reverse pass]
              + [alpha * W^T z,  z = (I - alpha M)^-T g]

    ⚠ DROPPING THE SECOND TERM WOULD NOT LOOK WRONG. Measured, the two are
    comparable — 140 against 749 at `l = 0` — so an implementation with
    only the first returns a plausible number rather than an obviously
    broken one. Hence a comparison against `m` independent forward driven
    solves, not a reasonableness check.

    ⚠ AND THE CIRCUIT HAS TO BE NONLINEAR FOR THIS TO MEAN ANYTHING. On a
    linear circuit the linearisation is time-INvariant and every `l != 0`
    is zero, so an `l = 1` comparison would be two noise-level numbers
    agreeing. The diode makes the sidebands real: measured
    `|H_1|/|H_0| = 0.50` and `|H_-2|/|H_0| = 0.13`.
    """
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e4)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    _pss, rows, fwds = _sideband_pair(cir, 1e-6, 1e6, 1, [l])

    assert np.abs(fwds[0]).max() > 1e-3, \
        'sideband l=%d is numerically absent (%.3g), so this comparison ' \
        'would be two zeros agreeing -- the circuit is not mixing' \
        % (l, np.abs(fwds[0]).max())
    rel = np.linalg.norm(rows[0] - fwds[0]) / np.linalg.norm(fwds[0])
    assert rel < 1e-11, \
        'sideband l=%d: the adjoint row disagrees with %d forward driven ' \
        'solves by %.3e' % (l, cir.n - 1, rel)


def test_sidebands_vanish_on_a_linear_circuit():
    """⚠ THE CHECK THAT CAUGHT A CONVENTION ERROR, and could only have.

    A linear circuit's linearisation is time-INvariant however hard it is
    driven, so it converts nothing: `H_l = 0` for every `l != 0`. Measured
    here at ~10 orders below `H_0`.

    The first implementation decomposed `y(t)` rather than
    `v(t) = y(t) exp(-j w t)`. That is self-consistent — it agreed with a
    forward reference written the same way to 1e-15 — but it smears a
    Dirichlet kernel across every `l` whenever `f` is not a multiple of
    `1/T`, and reported `|H_1|` comparable to `|H_0|` on this very circuit.
    Only a case whose answer is known independently separates the two.
    """
    cir = _adjoint_ladder(4)
    _pss, rows, _f = _sideband_pair(cir, 1e-3, 700.0, 1, [0, 1, -2, 3])
    h0 = np.abs(rows[0]).max()
    assert h0 > 1.0, 'H_0 is empty (%.3g); nothing is being measured' % h0
    for li, l in enumerate([0, 1, -2, 3][1:], start=1):
        rel = np.abs(rows[li]).max() / h0
        assert rel < 1e-8, \
            'a LINEAR circuit converted %.3g of H_0 into sideband l=%d. ' \
            'It has no time-varying linearisation to convert with, so this ' \
            'is the sideband convention, not the circuit' % (rel, l)


def test_H0_is_the_ac_transfer_function_and_converges_like_gear():
    """`H_0` against a reference the adjoint path cannot influence.

    On a linear circuit the `l = 0` row IS the LTI transfer function from
    each source to the output — an `(sC + G)` solve, no monodromy, no
    replay, no adjoint. That makes it the strongest available check of the
    whole `H_l` path, and it is the one the pnoise build should be able to
    lean on.

    ⚠ THE RATE IS THE ASSERTION. The residual at 60 points is 1.2e-03,
    which in isolation says nothing — it could be a defect of any size.
    Per doubling of the grid it falls 4.06x / 4.03x / 4.02x, which is
    Gear-2's own O(h^2) and therefore pure discretisation.
    """
    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric
    per, freq, k = 1e-3, 700.0, 1
    rels = []
    for npts in (60, 120, 240):
        cir = _adjoint_ladder(4)
        m = cir.n - 1
        _pss, rows, _f = _sideband_pair(cir, per, freq, k, [0], npts=npts,
                                        reltol=1e-12)
        irn = _pss.irefnode
        G = np.asarray(cir.G(np.zeros(cir.n)))
        Cm = np.asarray(cir.C(np.zeros(cir.n)))
        G, Cm = remove_row_col((G, Cm), irn, _pss.toolkit)
        sM = 2j * np.pi * freq * np.asarray(Cm) + np.asarray(G)
        ac = np.array([-np.linalg.solve(sM, np.eye(m)[:, i])[k]
                       for i in range(m)])
        rels.append(np.linalg.norm(rows[0] - ac) / np.linalg.norm(ac))

    for a, b in zip(rels, rels[1:]):
        ratio = a / b
        assert 3.4 < ratio < 4.6, \
            'H_0 against the AC transfer converges at %.2fx per doubling ' \
            '(%.4e -> %.4e), not the O(h^2) Gear-2 gives everywhere else. ' \
            'A wrong constant is invisible to a size check and obvious ' \
            'here' % (ratio, a, b)


def test_the_sideband_bound_is_hard():
    """⚠ THE TRUNCATION BOUND IS A REFUSAL, NOT A TOLERANCE.

    Okumura et al. eq. (32): the maximum frequency the analysis can speak
    about is the grid's own, so `|l| <= (w_max - w0)/ws`. Nothing aliases
    down from above what the grid represents. The abstract's "accumulated
    until their contributions become negligible" is a ratio test operating
    INSIDE that ceiling — an implementation carrying only the ratio test
    stops for the wrong reason, and on a coarse grid it would stop having
    summed harmonics the grid cannot represent.

    So this is refused rather than clamped: the remedy is a finer period
    grid, and silently returning the nearest representable sideband would
    answer a question the caller did not ask.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 40, maxiterations=40)
    fp = pss.factored_period()
    lmax = len(fp.steps) // 2
    pac = PAC(cir, toolkit=circuit.numeric)

    pac.adjoint_sideband_row(pss, 700.0, 1, lmax)      # at the edge: fine
    with pytest.raises(ValueError, match='Nyquist'):
        pac.adjoint_sideband_row(pss, 700.0, 1, lmax + 1)


def test_the_sideband_row_costs_one_solve_per_sideband():
    """The many-to-one property has to survive the distributed output.

    The injected reverse pass and the transposed solve both depend on `l`,
    so the cost is per SIDEBAND — but still not per SOURCE, which is the
    asymmetry pnoise is shaped by. This asserts the count, not the clock:
    this machine runs more than one agent.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(10)
    per, freq = 1e-3, 700.0
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 60, maxiterations=40)
    pac = PAC(cir, toolkit=circuit.numeric)
    rows, rows_info = pac.adjoint_sideband_row(pss, freq, 1, [0, 1, 2])
    assert rows.shape == (3, m)
    per_sideband = rows_info['matvecs'] / 3.0
    assert per_sideband < m, \
        'the sideband row took %.1f matvecs per sideband at m=%d; it is ' \
        'supposed to be independent of the number of sources' \
        % (per_sideband, m)


def test_the_pac_operator_is_singular_at_every_harmonic_of_an_oscillator():
    """⚠ AT DC *AND* AT EVERY HARMONIC, AND ONLY FOR AN OSCILLATOR.

    `I - exp(-j w T) M` is what PAC solves. At `w = k w0` the factor is 1
    and the operator is `I - M`, which an autonomous circuit's unit
    multiplier makes singular. Measured on van der Pol:

        offset/f0    0      0.25   0.5    0.75    1       2       3
        sigma_min  2.8e-11  0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

    and LINEAR in the distance to the nearest harmonic — 2.5e-1, 2.6e-2,
    2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way. The same sweep
    on a DRIVEN ladder never drops below 0.78.

    ⚠ IT IS PHYSICS, NOT CONDITIONING, which is why PAC refuses rather
    than tightening a tolerance: a perturbation at a harmonic is a
    perturbation ALONG the orbit, and an oscillator answers that with
    unbounded phase drift. There is no bounded periodic response to
    return, so a number there would be a wrong answer rather than an
    imprecise one.
    """
    circuit.default_toolkit = circuit.numeric

    _cir, pss = _solve_slow(None)                 # van der Pol, autonomous
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])

    def smin(ratio):
        a = np.exp(-2j * np.pi * ratio)
        return np.linalg.svd(np.eye(n) - a * M, compute_uv=False)[-1]

    for k in (0.0, 1.0, 2.0, 3.0):
        assert smin(k) < 1e-7, \
            'harmonic %g is not singular (%.3e); the unit multiplier has ' \
            'gone and so has the oscillator' % (k, smin(k))
    for r in (0.25, 0.5, 0.75):
        assert smin(r) > 0.1, \
            'offset %g is singular too (%.3e), so the singularity is not ' \
            'specific to the harmonics' % (r, smin(r))
    ## linear approach, not quadratic and not a cliff
    near = [smin(1.0 - d) for d in (1e-2, 1e-3, 1e-4)]
    for a, b in zip(near, near[1:]):
        assert 5 < a / b < 20, \
            'sigma_min falls %.1fx per decade of distance to the harmonic ' \
            '(%.3e -> %.3e), not the linear 10x' % (a / b, a, b)

    ## and PAC refuses there rather than returning something
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    f0 = 1.0 / pss.period
    with pytest.raises(ValueError, match='harmonic'):
        pac.adjoint_sideband_row(pss, f0, 0, 0)
    with quiet():
        pac.adjoint_sideband_row(pss, 0.37 * f0, 0, 0)   # off-harmonic: fine

    ## a DRIVEN circuit has no such structure and is not refused
    cir2 = _adjoint_ladder(3)
    pss2 = PSS(cir2, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    fp2 = pss2.factored_period()
    n2 = fp2.width
    M2 = np.column_stack([fp2.matvec(e) for e in np.eye(n2)])
    for k in (0.0, 1.0, 2.0):
        s2 = np.linalg.svd(np.eye(n2) - np.exp(-2j * np.pi * k) * M2,
                           compute_uv=False)[-1]
        assert s2 > 1e-3, \
            'a DRIVEN circuit is singular at harmonic %g (%.3e)' % (k, s2)
    PAC(cir2, toolkit=circuit.numeric).adjoint_sideband_row(
        pss2, 1.0 / pss2.period, 1, 0)             # not refused


def test_the_am_pm_split_is_exact_and_the_conjugate_is_load_bearing():
    """⚠ ONE CONJUGATE SEPARATES AM FROM PM, and `a ± b` looks just as right.

    The two sidebands COUNTER-ROTATE about the carrier phasor, so the sum
    traces an ellipse: the component along the carrier is amplitude
    modulation, perpendicular is phase modulation. Pure AM keeps the
    envelope on the carrier's axis, forcing `a = conj(b)`; pure PM keeps it
    perpendicular, `a = −conj(b)`. Hence `m_am = a + conj(b)` and
    `m_pm = a − conj(b)`, each vanishing exactly when the other case holds.

    ⚠ WITHOUT THE CONJUGATE the split still produces two numbers and they
    are wrong for any modulation whose sidebands are not real relative to
    the carrier — it reports a rotating ellipse as pure AM. This test pins
    that by checking the naive form does NOT vanish where the correct one
    does, so a "simplification" back to `a ± b` fails here rather than in
    somebody's phase-noise number.
    """
    rng = np.random.default_rng(0)
    for _ in range(5):
        a = complex(rng.standard_normal(), rng.standard_normal())

        m_am, m_pm = PAC.am_pm_indices(a, np.conj(a))        # pure AM
        assert abs(m_pm) < 1e-14 * max(abs(m_am), 1.0), \
            'pure AM leaked %.3e into the PM index' % abs(m_pm)
        m_am2, m_pm2 = PAC.am_pm_indices(a, -np.conj(a))     # pure PM
        assert abs(m_am2) < 1e-14 * max(abs(m_pm2), 1.0), \
            'pure PM leaked %.3e into the AM index' % abs(m_am2)

    ## a case where the naive split is unambiguously wrong: pure AM with a
    ## complex modulation phase. `a + b` is then NOT zero-PM.
    a = complex(0.3, 0.7)
    b = np.conj(a)                       # pure AM by construction
    _m_am, m_pm = PAC.am_pm_indices(a, b)
    assert abs(m_pm) < 1e-14
    assert abs(a - b) > 0.5, \
        'the naive `a - b` should be far from zero on this pure-AM case, ' \
        'which is exactly why the conjugate cannot be dropped'


def test_am_pm_on_a_driven_mixer_is_predominantly_am():
    """A diode detector converts amplitude to amplitude; it has no free phase.

    So the modulation a small signal imposes on the carrier should come out
    overwhelmingly AM — measured `|m_pm|/|m_am| = 0.005` — and should be
    stable across modulation offsets, since neither index has a `1/f`
    mechanism behind it here. (An oscillator is the opposite case: its
    phase response goes as `1/ω_m`, so PM dominates near the carrier.)
    """
    circuit.default_toolkit = circuit.numeric
    cir = _diode_mixer()
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-6, timestep=1e-6 / 120, maxiterations=40)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index(2)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / pss.period

    ## the carrier has to exist before it can be modulated
    C1 = pac.carrier_phasor(pss, k, 1)
    assert abs(C1) > 1e-3, 'no fundamental at the output (%s)' % C1

    ratios = []
    for r in (0.3, 0.1, 0.03):
        with quiet():
            am, pm = pac.am_pm(pss, r * f0, k, 1, sweeptype='relative')
        i = int(np.argmax(np.abs(am)))
        assert abs(am[i]) > 1.0, \
            'the AM index is %.3e -- nothing is being modulated' % abs(am[i])
        ratios.append(abs(pm[i]) / abs(am[i]))

    assert max(ratios) < 0.05, \
        'a diode detector came out %.4f PM/AM; it has no free phase to ' \
        'modulate' % max(ratios)
    assert max(ratios) / min(ratios) < 1.2, \
        'the PM/AM ratio moved %.2fx across offsets (%s); neither index has ' \
        'a 1/f mechanism here, so it should be flat' \
        % (max(ratios) / min(ratios), ratios)


def test_am_pm_refuses_a_harmonic_with_no_carrier():
    """AM/PM of nothing is not a small number, it is undefined."""
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    pac = PAC(cir, toolkit=circuit.numeric)
    with pytest.raises(ValueError, match='no component at harmonic'):
        pac.am_pm(pss, 137.0, 1, harmonic=57,
                  sweeptype='relative')     # far above anything present


def test_the_deflated_solve_removes_the_harmonic_singularity():
    """⚠ A7 — THE POLE IS THE ANSWER'S, AND IT IS NOW CARRIED ANALYTICALLY.

    At `α = 1` the operator `I − αM` is singular by the unit multiplier,
    and the solution really does diverge: an oscillator's phase response
    goes as `1/Δf`. What is wrong is *computing* a `1/ε` quantity through a
    system whose conditioning is also `1/ε` — the answer is genuinely large
    and the digits are genuinely gone.

    Bordering with the tangent and the PPV — both of which `ppv()` already
    returns — makes the border variable bounded, `s = (vᵀb)/(vᵀu)`, and the
    pole comes back through `u/(1 − α)` in closed form.

    Measured here: the plain operator's `σ_min` tracks the offset over nine
    decades while the bordered one is FLAT, and the solutions agree to
    ~1e-11 where the plain solve is still trustworthy.

    ⚠ THE TEST ASSERTS BOTH HALVES. Flat conditioning alone would be
    satisfied by an operator that had stopped solving the right problem, so
    agreement with the plain solve — in the regime where the plain solve
    can be believed — is what says it is still the same equation.
    """
    _cir, pss = _solve_slow(None)                # van der Pol, autonomous
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    rng = np.random.default_rng(0)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)

    smins, rels = [], []
    for r in (1e-2, 1e-4, 1e-6):
        alpha = np.exp(-2j * np.pi * r)
        with quiet():
            y = pac._deflated_solve(pss, alpha, b)
        smins.append(np.linalg.svd(np.eye(n) - alpha * M,
                                   compute_uv=False)[-1])
        yp = np.linalg.solve(np.eye(n) - alpha * M, b)
        rels.append(np.linalg.norm(y - yp) / np.linalg.norm(y))

    ## the plain operator really is degrading, or this proves nothing
    assert smins[0] / smins[-1] > 1e3, \
        'the plain operator only degraded %.1fx over four decades of ' \
        'offset (%s); the singularity being removed is not being exercised' \
        % (smins[0] / smins[-1], smins)

    ## and where the plain solve can still be believed, the two agree
    assert rels[0] < 1e-7, \
        'the deflated solve disagrees with the plain one by %.3e at an ' \
        'offset where the plain one is well conditioned — it is solving a ' \
        'different equation' % rels[0]

    ## the physical pole survives: |y| grows as 1/df
    mags = []
    for r in (1e-3, 1e-4, 1e-5):
        with quiet():
            mags.append(np.linalg.norm(
                pac._deflated_solve(pss, np.exp(-2j * np.pi * r), b)))
    for a, b_ in zip(mags, mags[1:]):
        assert 8.0 < b_ / a < 12.0, \
            'the response grows %.2fx per decade closer to the harmonic, ' \
            'not the 10x of a 1/df pole — the pole has been removed from ' \
            'the ANSWER rather than from the conditioning' % (b_ / a)


def test_the_deflated_solve_still_refuses_an_exact_harmonic():
    """Because the physical response there is unbounded.

    `1/(1 − α)` is a division by zero at a harmonic. The pole is removed
    from the CONDITIONING, not from the answer, so an exact harmonic is
    still a question with no finite answer.
    """
    _cir, pss = _solve_slow(None)
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    n = pss.factored_period().width
    with pytest.raises(ValueError, match='EXACT harmonic'):
        pac._deflated_solve(pss, 1.0 + 0.0j, np.ones(n))


def test_the_deflated_solve_works_transposed_too():
    """The adjoint rows need `(I − αMᵀ)`, whose borders swap.

    Its null space is spanned by the PPV and its left null space by the
    tangent, so `u` and `v` exchange roles. Checked against a dense
    transposed solve where that is still trustworthy.
    """
    _cir, pss = _solve_slow(None)
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    rng = np.random.default_rng(3)
    d = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    alpha = np.exp(-2j * np.pi * 1e-2)
    with quiet():
        x = pac._deflated_solve(pss, alpha, d, transposed=True)
    xp = np.linalg.solve(np.eye(n) - alpha * M.T, d)
    rel = np.linalg.norm(x - xp) / np.linalg.norm(x)
    assert rel < 1e-7, \
        'the transposed deflated solve disagrees with a dense transposed ' \
        'solve by %.3e; the borders are probably not swapped' % rel


def test_the_sideband_response_matches_the_analytic_cos2_coefficients():
    """⚠⚠ PAC'S ANSWER IS INDEXED BY SIDEBAND, and every entry is known here.

    Kundert: *"for a single output frequency there may be many transfer
    functions from a single input"*. `SidebandResponse` carries all of
    them; a caller who takes one coefficient and calls it "the gain" has
    silently picked one.

    On the `cos²` mixer every entry is analytic — the transconductance's
    Fourier coefficients are `H₀ = ½`, `H_±2 = ¼`, and **`H_±1 = 0`**:

        l    f_in (Hz)   |H|/gain    expected
        −2      4500.0   0.25000000    0.25
        −1      3500.0   0.00000000    0.00
        +0      2500.0   0.50000000    0.50
        +1      1500.0   0.00000000    0.00
        +2       500.0   0.25000000    0.25

    ⚠ THE ZEROS ARE THE DISCRIMINATING PART. A sideband index off by one
    would move weight onto `l = ±1`, and the magnitudes alone would still
    look like a plausible mixer. `rejection_db(0, 2)` comes out at
    **6.0206 dB = 20·log₁₀(2)** exactly, and the `±2` pair is symmetric to
    1e-15.

    ⚠⚠ AND EACH SIDEBAND IS FED FROM ITS OWN INPUT BAND — `f_in = f_out −
    l·f₀`, five different frequencies for one output. That is why image
    rejection is not `H_l` against `H_−l` at a single input: computing it
    that way gives a plausible number for a different quantity.
    """
    _cir, pss, pac, d, gain2 = _cos2_mixer()
    cir = _cir
    g = np.sqrt(gain2)
    irn = pss.irefnode
    src = cir.get_node_index('vin')
    src = src - 1 if src > irn else src

    r = pac.mixer_response(pss, 2500.0, d, sidebands=(-2, -1, 0, 1, 2))
    assert r.sidebands == [-2, -1, 0, 1, 2]
    f0 = 1.0 / float(pss.period)
    for l in r.sidebands:
        assert abs(r.input_frequency(l) - (2500.0 - l * f0)) < 1e-9, \
            'sideband %d is fed from %.6g Hz, not f_out - l*f0' \
            % (l, r.input_frequency(l))

    expect = {0: 0.5, 2: 0.25, -2: 0.25, 1: 0.0, -1: 0.0}
    for l, want in expect.items():
        got = abs(r.transfer(l, src)) / g
        if want == 0.0:
            assert got < 1e-9, \
                'sideband %d should be a NULL for a cos^2 transconductance ' \
                'and reads %.3e; weight there means the sideband index is ' \
                'off' % (l, got)
        else:
            assert abs(got - want) < 1e-6, \
                'sideband %d reads %.8f against the analytic %.4f' \
                % (l, got, want)

    assert abs(r.rejection_db(0, 2, src) - 20.0 * np.log10(2.0)) < 1e-6
    assert abs(r.rejection_db(2, -2, src)) < 1e-9, \
        'the +-2 pair is not symmetric (%.3e dB)' % r.rejection_db(2, -2, src)
    assert r.rejection_db(0, 1, src) > 100.0, \
        'rejection against a null should be enormous, not %.3f' \
        % r.rejection_db(0, 1, src)


def test_the_sideband_row_is_indexed_by_SOURCE_not_by_the_output_direction():
    """⚠⚠ THE MISLABELLING THIS OBJECT EXISTS TO PREVENT, COMMITTED WHILE
    BUILDING IT — so it is pinned rather than merely remembered.

    `adjoint_sideband_row` returns a row over **sources**. Probing it at
    the index where the OUTPUT direction peaks returns the direct path
    from a current injected at the output node — a correct answer to a
    different question, and on this circuit a very convincing one:

        source index 1 (`vin`, through the mixer):  l=0 → ½·gain, l=±2 → ¼·gain
        source index 2 (`out`, direct through Rout): l=0 → 1·gain, l=±2 → ~0

    Read at index 2 the mixer looks like it has **no** sideband response
    at all and unity conversion gain. Every number is right; the column is
    wrong. Nothing in the shapes disagrees, because both are length-`m`
    rows of plausible magnitudes.

    That is why `SidebandResponse.transfer` takes the source explicitly
    and has no default: there is no sensible one, and the convenient
    guess is the output direction.
    """
    cir, pss, pac, d, gain2 = _cos2_mixer()
    g = np.sqrt(gain2)
    irn = pss.irefnode
    i_in = cir.get_node_index('vin')
    i_in = i_in - 1 if i_in > irn else i_in
    i_out = cir.get_node_index('out')
    i_out = i_out - 1 if i_out > irn else i_out
    assert i_in != i_out
    assert int(np.argmax(np.abs(d))) == i_out, \
        'the output direction no longer peaks at the output node, so the ' \
        'confusion this test documents is not reproducible'

    rows = np.asarray(pac.adjoint_sideband_row(pss, 500.0, d,
                                               sidebands=(0, 2))[0])
    thru = np.abs(rows[:, i_in]) / g
    direct = np.abs(rows[:, i_out]) / g
    assert abs(thru[0] - 0.5) < 1e-6 and abs(thru[1] - 0.25) < 1e-6, \
        'the mixer path no longer reads 1/2 and 1/4 (%s)' % thru
    assert abs(direct[0] - 1.0) < 1e-6 and direct[1] < 1e-9, \
        'the direct output path no longer reads 1 and 0 (%s); the two ' \
        'columns must stay distinguishable or this test proves nothing' \
        % direct


def test_mixer_response_refuses_a_negative_input_band():
    """A sideband whose input band would be negative is refused, not folded.

    For `f_out < l·f₀` the input frequency comes out negative. That band
    is the conjugate of `|f_in|`, so taking the absolute value quietly
    would return a right magnitude under a wrong label — the same class of
    error as reading the wrong column, and equally invisible.
    """
    _cir, pss, pac, d, _g = _cos2_mixer()
    ok = pac.mixer_response(pss, 2500.0, d, sidebands=(1,))
    assert ok.input_frequency(1) > 0
    with pytest.raises(ValueError, match='which is negative'):
        pac.mixer_response(pss, 300.0, d, sidebands=(1,))


def test_the_in_tree_gmres_matches_scipy_and_keeps_its_hessenberg():
    """⚠ WRITTEN RATHER THAN IMPORTED, AND SPEED IS NOT THE REASON.

    `_arnoldi_gmres` exists for two things SciPy's cannot do:

    **It keeps `H`.** García, Romero & Acha read the Floquet multipliers
    off exactly this matrix — Ritz values `θ` of `I − M` map back as
    `λ = 1 − θ` — so a GMRES that discards its Hessenberg matrix throws
    away the spectrum it just computed.

    **It judges by residual, not by a status flag.** SciPy returns
    `info = 4` on a HAPPY breakdown — the Krylov space exhausted *because
    the answer is exact* — which turned an exact AM/PM result into a
    `RuntimeError` and is why `PAC._gmres_checked` had to overrule it.
    Here the breakdown is detected where it happens.

    Gated three ways: agreement with SciPy to ~1e-15 on well-conditioned
    systems; Ritz values of `H` recovering a known spectrum to 1e-15 at
    `k = n`; and reorthogonalisation earning its cost — **5 orders** on
    the Ritz values (3.51e-10 → 4.48e-15 at a spectral spread of 1e2),
    because modified Gram-Schmidt loses orthogonality as the basis grows
    and a drifted basis gives multipliers **the residual cannot see are
    wrong**.
    """
    import scipy.sparse.linalg as spla
    from pycircuit.circuit.shooting._numerics import _arnoldi_gmres
    rng = np.random.default_rng(0)

    for n in (4, 12, 40):
        A = rng.standard_normal((n, n)) + n * np.eye(n)
        b = rng.standard_normal(n)
        x1, relres, _H, _k = _arnoldi_gmres(lambda v: A @ v, b, rtol=1e-13)
        x2, info = spla.gmres(A, b, rtol=1e-13, restart=min(n, 50),
                              maxiter=min(n, 200))
        assert info == 0, 'scipy failed on the reference system'
        assert np.linalg.norm(x1 - x2) / np.linalg.norm(x1) < 1e-10, \
            'n=%d: the in-tree solve disagrees with scipy' % n
        assert (np.linalg.norm(A @ x1 - b) / np.linalg.norm(b)
                < max(1e3 * 1e-13, 1e-10)), \
            'n=%d: relres reported %.3e but the true residual disagrees' \
            % (n, relres)

    ## H carries the spectrum
    n = 10
    ev = np.linspace(0.1, 0.9, n)
    Qr, _r = np.linalg.qr(rng.standard_normal((n, n)))
    A = Qr @ np.diag(ev) @ Qr.T
    b = rng.standard_normal(n)
    _x, _rr, H, k = _arnoldi_gmres(lambda v: A @ v, b, rtol=1e-15)
    assert k == n, 'the full basis was not built (k = %d)' % k
    ritz = np.sort(np.real(np.linalg.eigvals(H)))
    assert np.max(np.abs(ritz - ev)) < 1e-12, \
        'the Ritz values of H do not recover the spectrum (max err %.3e); ' \
        'H is the whole reason this function exists' \
        % float(np.max(np.abs(ritz - ev)))

    ## reorthogonalisation is not decoration
    ev2 = np.logspace(0, 2, 30)
    Q2, _r2 = np.linalg.qr(np.random.default_rng(7).standard_normal((30, 30)))
    A2 = Q2 @ np.diag(ev2) @ Q2.T
    b2 = np.random.default_rng(7).standard_normal(30)
    errs = {}
    for ro in (False, True):
        _x2, _r3, H2, k2 = _arnoldi_gmres(lambda v: A2 @ v, b2,
                                          rtol=1e-15, reortho=ro)
        assert k2 == 30
        errs[ro] = float(np.max(np.abs(
            np.sort(np.real(np.linalg.eigvals(H2))) - ev2) / ev2))
    assert errs[True] < errs[False] / 100.0, \
        'reorthogonalisation buys only %.1fx on the Ritz values (%.3e -> ' \
        '%.3e); if it has stopped mattering the extra pass should go, and ' \
        'if it has stopped WORKING the multipliers are wrong in a way the ' \
        'residual cannot see' % (errs[False] / errs[True],
                                 errs[False], errs[True])


def _ampm_mixer():
    """A NONLINEAR driven mixer with three non-reference nodes.

    ⚠ A LINEAR CIRCUIT IS USELESS HERE and the first attempt used one. An
    LTI circuit does no mixing, so the carrier has no sidebands and both
    modulation indices come out at ~1e-12 — noise, against which any
    invariance claim is noise-over-noise. The nonlinearity is what makes
    the operating point genuinely time-varying.
    """
    cir = SubCircuit()
    for nn in ('a', 'b', 'c'):
        cir.add_node(nn)
    cir['vs'] = VSin('a', gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    cir['R1'] = R('a', 'b', r=1e4)
    cir['D'] = Diode('b', gnd)
    cir['C1'] = C('b', gnd, c=1e-12)
    cir['R2'] = R('b', 'c', r=1e4)
    cir['C2'] = C('c', gnd, c=1e-12)
    return cir


def _ampm_at(refname, fm=5e4, npts=200):
    circuit.default_toolkit = circuit.numeric
    cir = _ampm_mixer()
    rn = gnd if refname == 'gnd' else cir.get_node(refname)
    pss = PSS(cir, method='gear', reltol=1e-12, irefnode=rn)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-6, timestep=1e-6 / npts, refnode=rn,
                  maxiterations=60)
    assert pss.converged, refname
    irn = pss.irefnode
    m = cir.n - 1

    def red(nm):
        k = cir.get_node_index(nm)
        return None if k == irn else (k - 1 if k > irn else k)

    def vec(p, q):
        v = np.zeros(m)
        for nm, sg in ((p, 1.0), (q, -1.0)):
            i = red(nm)
            if i is not None:
                v[i] += sg
        return v

    out = vec('b', 'c')
    src = vec('c', 'b')
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        am, pm = pac.am_pm(pss, fm, out, harmonic=1, sweeptype='relative')
        car = pac.carrier_phasor(pss, out, 1)
    am = np.asarray(am).ravel()
    pm = np.asarray(pm).ravel()
    return complex(np.sum(am * src)), complex(np.sum(pm * src)), complex(car)


def test_am_pm_accepts_a_direction_not_only_a_node_index():
    """⚠ `am_pm` COULD NOT EXPRESS A DIFFERENTIAL OUTPUT AT ALL, and that
    was an API asymmetry with a real consequence.

    `pnoise`, `adjoint_transfer_row` and `adjoint_sideband_row` all take a
    direction vector `d`. `am_pm` and `carrier_phasor` did `int(output)`
    and took a node index. So a differential observable — `v(b) − v(c)` —
    was not expressible, and for an oscillator the output of interest is
    very often differential.

    Fixed by `_output_waveform_row`, which accepts either. The integer
    form still works, so callers naming a node are unaffected; this pins
    both, and pins that they agree where they overlap.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _ampm_mixer()
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-6, timestep=1e-6 / 200, refnode=gnd,
                  maxiterations=60)
    assert pss.converged
    irn = pss.irefnode
    kb = cir.get_node_index('b')
    kb = kb - 1 if kb > irn else kb
    pac = PAC(cir, toolkit=circuit.numeric)

    ## the integer form still works and agrees with the equivalent vector
    e_b = np.zeros(cir.n - 1)
    e_b[kb] = 1.0
    with quiet():
        c_int = pac.carrier_phasor(pss, kb, 1)
        c_vec = pac.carrier_phasor(pss, e_b, 1)
    assert abs(c_int - c_vec) <= 1e-12 * abs(c_int), \
        'the index and single-entry-vector forms disagree (%r vs %r)' \
        % (c_int, c_vec)

    ## and a differential output is now expressible AND different
    kc = cir.get_node_index('c')
    kc = kc - 1 if kc > irn else kc
    diff = np.zeros(cir.n - 1)
    diff[kb], diff[kc] = 1.0, -1.0
    with quiet():
        c_diff = pac.carrier_phasor(pss, diff, 1)
        am, pm = pac.am_pm(pss, 5e4, diff, harmonic=1, sweeptype='relative')
    assert abs(c_diff - c_int) > 1e-3 * abs(c_int), \
        'the differential carrier equals the single-ended one, so this ' \
        'fixture does not exercise the new path'
    assert np.asarray(am).size == cir.n - 1
    assert np.all(np.isfinite(np.asarray(am)))


def test_am_pm_is_invariant_under_a_change_of_reference_node():
    """⚠⚠ KÄRTNER §3.2 — the property an AM/PM split must have, and what
    this test does and does NOT establish.

    Under a linear change of state variables the Floquet basis transforms
    as `u' = Au`, `v'ᵀ = vᵀA⁻¹`, and Kärtner obtains *"the same equation
    for the time shift θ(t) … therefore the separation in amplitude and
    phase is **independent of the co-ordinate system used** … this is by
    no means a trivial result, since there are **arbitrarily many other
    definitions of amplitude and phase which seem to be more illustrative
    but do not have this invariance**, and therefore a change of
    co-ordinates also **transforms a part of phase noise into amplitude
    noise**."*

    A change of reference node is such a transformation — a genuine shear,
    not a permutation — and the split of a *differential* observable comes
    out invariant. Measured across `gnd`, `a` and `c`:

        |C₁| identical to 4e-16,  m_am to 1.5e-12,  m_pm to 1.2e-10

    with `m_am ≈ 28.3` and `m_pm ≈ 1.45`, so these are real numbers rather
    than agreement between two zeros.

    ⚠⚠ WHAT IT DOES NOT ESTABLISH, STATED BECAUSE THE CITATION INVITES THE
    STRONGER READING. The sidebands of a *fixed physical observable* are
    themselves coordinate-free, so **any** function of them would pass
    this — including the plausible-but-wrong definitions Kärtner warns
    about. What this pins is that the IMPLEMENTATION handles the reference
    row correctly: `_output_waveform_row`'s reinsertion, the reduced-index
    mapping, and the direction contraction. That is worth holding and it
    is not Kärtner's theorem.

    Discriminating among AM/PM *definitions* needs a different experiment
    — one that transforms the state basis rather than the observable's
    representation — and is not built.
    """
    res = {}
    for rn in ('gnd', 'a', 'c'):
        res[rn] = _ampm_at(rn)

    base = res['gnd']
    assert abs(base[0]) > 1.0 and abs(base[1]) > 0.1, \
        'the modulation indices are ~0 (%r), so this fixture has stopped ' \
        'being nonlinear and the invariance below would be noise against ' \
        'noise' % (base,)
    for rn in ('a', 'c'):
        for j, nm, tol in ((0, 'm_am', 1e-9), (1, 'm_pm', 1e-7),
                           (2, 'C1', 1e-12)):
            rel = abs(res[rn][j] - base[j]) / abs(base[j])
            assert rel < tol, \
                'refnode %s moved %s by %.3e relative; the split of a ' \
                'differential observable must not depend on which row was ' \
                'eliminated' % (rn, nm, rel)


def test_the_pac_adjoint_surfaces_run_under_every_integrator():
    """⚠ The other two refusals B8 left standing: `adjoint_transfer_row`
    and `adjoint_sideband_row` both keyed on `fp.kind` and both said the
    transposed replay was "implemented for the solved-history map only" —
    a sentence that stopped being true when the plain recursion shipped.

    Both consume the reverse pass through `collect=True` and `inject=`,
    so the plain replay had to grow those too; `_forced_replay_transposed`
    reads `ts[j]`, which is the transposed solve at step `j` under BOTH
    recursions.

    Measured against gear on a driven RLC (relative, at 200 then 400
    points):

        trap   adjoint 7.5e-05 -> 3.8e-05    sideband 2.1e-06 -> 1.2e-06
        euler  adjoint 1.7e-04 -> 8.6e-05    sideband 1.7e-04 -> 8.6e-05

    ⚠ `trap` converging at `O(h)` rather than `O(h²)` is NOT a defect here
    and not a surprise: the manufacturing step costs PAC an order already,
    which `test_pac_order_is_lost_to_the_manufacturing_step` pins
    independently with `x0_unknown` as the switch.
    """
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    Lv, Cv, Rs = 1e-3, 1e-9, 10.0
    per = 2.0 * np.pi * np.sqrt(Lv * Cv)

    def build():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['r'] = R('a', 'b', r=Rs)
        c['l'] = L('b', gnd, L=Lv)
        c['c1'] = C('b', gnd, c=Cv)
        c['n'] = IS('b', gnd, i=0.0, noisePSD=1e-18)
        return c

    freq = 0.31 / per
    got = {}
    for method in ('gear', 'trap', 'euler'):
        for npts in (200, 400):
            cir = build()
            pss = PSS(cir, method=method, reltol=1e-11)
            with quiet(AccuracyWarning):
                pss.solve(period=per, timestep=per / npts,
                          x0=np.zeros(cir.n - 1), maxiterations=100,
                          x0_unknown=False)
            assert pss.converged
            ob = [str(nd) for nd in cir.nodes].index('b')
            pac = PAC(cir)
            with quiet():
                ar = pac.adjoint_transfer_row(pss, freq, ob)[0]
                sr = pac.adjoint_sideband_row(pss, freq, ob, sidebands=0)[0]
            got[(method, npts)] = (float(np.linalg.norm(ar)),
                                   float(np.linalg.norm(sr)))

    for npts in (200, 400):
        g = got[('gear', npts)]
        for method in ('trap', 'euler'):
            r = got[(method, npts)]
            for k, name in ((0, 'adjoint'), (1, 'sideband')):
                rel = abs(r[k] - g[k]) / abs(g[k])
                assert rel < 1e-3, \
                    '%s/%d: the %s row disagrees with gear by %.3e' \
                    % (method, npts, name, rel)

    ## and it must IMPROVE with refinement, or the agreement above is
    ## accidental rather than convergent
    for method in ('trap', 'euler'):
        a = abs(got[(method, 200)][0] - got[('gear', 200)][0])
        b = abs(got[(method, 400)][0] - got[('gear', 400)][0])
        assert b < a, \
            '%s: the adjoint row does not converge toward gear ' \
            '(%.3e then %.3e)' % (method, a, b)


def test_pac_sweep_recycling_makes_matvecs_INDEPENDENT_of_sweep_length():
    """⚠ THE RECYCLING IS BUILT AND DEFAULT-ON; THIS PINS WHAT IT BUYS.

    `_solve_subspace` shares one Krylov basis across the whole sweep, on
    Telichevesky's Theorem 1: `A(alpha) = I - alpha M`, so
    `span{r, Ar, A^2 r, ...} = span{r, Mr, M^2 r, ...}` for EVERY alpha.
    The basis is frequency-independent; each frequency then costs a small
    dense least-squares over it.

    Measured against `recycle=False` on a driven RLC with an RC ladder:

        m    K    mv recycle   mv each   ratio    t recyc   t each
        4    4         5          12      2.4x     0.149     0.106
        4   64         5         192     38.4x     0.787     1.702
       18    4         8          24      3.0x     0.270     0.179
       18   64        10         384     38.4x     0.974     2.862

    ⚠⚠ MATVECS ARE FLAT IN `K` AND THE WALL CLOCK IS NOT -- 38x against
    2.9x. That gap is the useful part: the Krylov solve has already been
    removed from the sweep's cost, and what remains is the ONE FORCED
    REPLAY PER FREQUENCY outside it, which recycling cannot touch. Anyone
    optimising this sweep further should go after the replays, not the
    linear solve.

    ⚠ AND IT IS A LOSS ON SHORT SWEEPS: at `K = 4` recycling costs MORE
    wall clock than solving each (0.149 against 0.106) despite using
    fewer matvecs, because the dense least-squares over the shared basis
    dominates. The win needs roughly `K >= 8`. That is not an argument
    against the default -- a 4-point PAC sweep is not where time goes --
    but it is why the ratio must be read on matvecs and length, not on a
    single timing.
    """
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    Lv, Cv, Rs = 1e-3, 1e-9, 10.0
    per = 2.0 * np.pi * np.sqrt(Lv * Cv)

    def build():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['r'] = R('a', 'b', r=Rs)
        c['l'] = L('b', gnd, L=Lv)
        c['c1'] = C('b', gnd, c=Cv)
        return c

    cir = build()
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pss.solve(period=per, timestep=per / 200,
                  x0=np.zeros(cir.n - 1), maxiterations=100,
                  x0_unknown=False)
    assert pss.converged

    seen = {}
    for K in (4, 32):
        freqs = np.linspace(0.05, 0.45, K) / per
        out = {}
        for flag in (True, False):
            pac = PAC(cir)
            with quiet():
                res = pac.solve(pss, freqs, recycle=flag)
            out[flag] = (res.info['matvecs'], np.asarray(res.x))
        ## the two routes must agree -- recycling minimises the TRUE
        ## residual over the shared span, so it is never the worse answer
        a, b = out[True][1], out[False][1]
        scale = float(np.max(np.abs(b)))
        rel = float(np.max(np.abs(a - b))) / max(scale, 1e-300)
        assert rel < 1e-9, \
            'K=%d: recycling changed the answer by %.3e' % (K, rel)
        seen[K] = (out[True][0], out[False][0])

    ## ⚠ THE UNRECYCLED COUNT MUST SCALE WITH THE SWEEP AND THE RECYCLED
    ## ONE MUST NOT. That contrast is the claim; an absolute count is a
    ## machine and tolerance detail.
    r4, e4 = seen[4]
    r32, e32 = seen[32]
    assert e32 > 6 * e4, \
        'solving each frequency no longer scales with sweep length ' \
        '(%d at K=4, %d at K=32)' % (e4, e32)
    assert r32 < 3 * r4, \
        'the recycled basis is no longer shared across the sweep -- its ' \
        'matvec count grew from %d to %d when the sweep grew 8x, which ' \
        'means the frequency-independence of the Krylov space has been ' \
        'lost' % (r4, r32)
    assert e32 > 5 * r32, \
        'recycling no longer saves matvecs at K=32 (%d against %d)' \
        % (e32, r32)


def test_pac_reports_sidebands_at_the_right_frequencies_and_conjugates_the_fold():
    """✅ TWO REPORTING DEFECTS IN `PAC.solve`'S TAIL, found by an external
    reference-simulator cross-check (2026-09-05), neither in the solve.  (a) The DFT was
    taken over `fp.times`, `[0, T]` INCLUSIVE, so the last sample repeated
    the first and `dt = T/(N-1)` put the sidebands at `f0 (N-1)/N`:
    109 500 / 89 500 Hz for 110 000 / 90 000 at N = 200 (measured before
    the fix: 109 500 / 109 750 / 109 875 at 200 / 400 / 800).  It cost an
    order, O(h) for O(h^2).  (b) `|sb + f|` folded a negative sideband
    frequency to positive and left the coefficient alone; the physical
    response there is the CONJUGATE.  Both were invisible on every
    earlier PAC gate, whose `v(t)` is constant over the period (a
    constant's DFT is exact for any window): it takes a CONVERTING
    circuit -- the switched capacitor -- to see them.  The proof is the
    adjoint: `adjoint_sideband_row` never calls `freq_analysis`, and after
    the fix the reported coefficients equal `H_l . u_ac` to 1e-15 at
    every grid, with `l = -1` the conjugate of `H_{-1} . u_ac`.
    """
    circuit.default_toolkit = circuit.numeric
    fclk, fin, T = 100e3, 10e3, 1e-5
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    cir.add_node('ck')
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0, vac=1.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                           vth=0.0, vs=50e-3)
    cir['C0'] = C('out', gnd, c=100e-12)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / 200, x0=np.zeros(cir.n - 1),
                  maxiterations=100)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    res = pac.solve(pss, freqs=[fin])
    fout = np.asarray(res.sweep_values, dtype=float)
    X = np.asarray(res.x)
    io = [str(n) for n in cir.nodes].index('out')
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), pss.irefnode,
                             circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    H = np.asarray(pac.adjoint_sideband_row(pss, fin, io,
                                            sidebands=[0, 1, -1])[0])
    for li, l in enumerate((0, 1, -1)):
        f_phys = abs(fin + l * fclk)
        k = int(np.argmin(np.abs(fout - f_phys)))
        assert abs(fout[k] - f_phys) < 1e-6 * fclk, \
            'sideband l=%d reported at %.1f Hz, not %.1f' % (l, fout[k], f_phys)
        x = complex(X[io, k])
        h = complex(H[li] @ u_ac)
        if fin + l * fclk < 0:
            h = np.conj(h)
        assert abs(x - h) < 1e-12 * abs(h), \
            'sideband l=%d: reported %r against the adjoint %r' % (l, x, h)


def test_the_sideband_gate_rejects_the_endpoint_and_the_unconjugated_fold():
    """MUTATION CHECK (P3): the two `PAC.solve` reporting defects, injected,
    and the gate's own criterion (agreement with `adjoint_sideband_row`)
    shown to reject each.  Both were invisible on a circuit whose `v(t)`
    is constant; this uses a converting circuit so the mutations bite.

    (a) ENDPOINT: taking the DFT over the `[0, T]`-inclusive grid puts the
    sidebands at `f0 (N-1)/N`, not `f0`, so a frequency-placement gate
    catches it.  (b) CONJUGATE: reporting a negative-fold coefficient
    without conjugating disagrees with `H_l . u_ac`, while the conjugate
    agrees -- the equality gate rejects the mutation and accepts the fix.
    """
    from pycircuit.circuit.shooting import freq_analysis
    circuit.default_toolkit = circuit.numeric
    fclk, fin, T = 100e3, 10e3, 1e-5
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    cir.add_node('ck')
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0, vac=1.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                           vth=0.0, vs=50e-3)
    cir['C0'] = C('out', gnd, c=100e-12)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / 200, x0=np.zeros(cir.n - 1),
                  maxiterations=100)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    io = [str(n) for n in cir.nodes].index('out')
    fp = pss.factored_period()
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), pss.irefnode,
                             circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    res = pac.solve(pss, freqs=[fin])
    fout = np.asarray(res.sweep_values, dtype=float)
    H = np.asarray(pac.adjoint_sideband_row(pss, fin, io, sidebands=[1, -1])[0])
    ## the FIX places l=+1 at 110 kHz and l=-1 at 90 kHz, and l=-1 equals
    ## conj(H_{-1}.u) -- the gate.  Confirm the gate is met, then that each
    ## mutation would break it.
    X = np.asarray(res.x)
    kpos = int(np.argmin(np.abs(fout - (fin + fclk))))
    kneg = int(np.argmin(np.abs(fout - abs(fin - fclk))))
    assert abs(fout[kpos] - (fin + fclk)) < 1e-6 * fclk        # placement OK
    assert abs(complex(X[io, kneg]) - np.conj(complex(H[1] @ u_ac))) \
        < 1e-9 * abs(complex(H[1] @ u_ac))                      # conjugate OK
    ## (a) endpoint mutation: the inclusive-window spacing is f0 (N-1)/N
    N = len(fp.steps)
    bad_spacing = fclk * (N - 1) / N
    assert abs(bad_spacing - fclk) > 1e-4 * fclk, \
        'the endpoint mutation must shift the sideband spacing; N=%d' % N
    ## (b) conjugate mutation: the UN-conjugated coefficient disagrees
    unconj = complex(X[io, kneg])            # the reported (fixed) value
    ## the mutation would report conj of the correct one; show they differ
    assert abs(unconj - np.conj(unconj)) > 1e-3 * abs(unconj), \
        'l=-1 has a real-only coefficient here, so the conjugate mutation ' \
        'is invisible on this fixture -- pick one with a phase'


def test_the_oscillator_am_pm_rows_are_the_isf_dc_term_and_vanish_by_half_wave_symmetry():
    """2026-09-08: the "~1e-12 sideband rows on an oscillator, cause not
    established" caveat, established.  `am_pm(pss, freq, output)` is the
    p = 0 band of the noise split: a source at BASEBAND `freq` reaching
    the carrier sideband.  A baseband current moves the oscillator's PHASE
    through the PPV's DC coefficient (Hajimiri-Lee's c_0, the 1/f^3
    up-conversion term), and a half-wave-symmetric orbit -- odd
    nonlinearity, `u(t + T/2) = -u(t)` -- has none: the same symmetry zero
    the coloured-upconversion gate records for Gamma.  Measured on
    `_lc_osc`, |m_pm| at 1e-3 f0:  a = 0: 1.24e-08;  a = 0.05: 46.6;
    a = 0.25: 231 -- a jump of 4e9 on breaking the symmetry, then LINEAR
    in `a` (470 per unit at both).  The direct p = 1 rows meanwhile agree
    with pnoise at every offset (1/sqrt 2 of sqrt(S/psd) near the carrier,
    where the image band carries the other half; 1.000 far out), and the
    three-leg chain puts the split on the certified Lorentzian's absolute
    scale, so nothing is small for an unestablished reason: the rows are
    right, and the symmetric fixture measures a zero.  ⚠ My estimate of
    the lifted magnitude from a tank-impedance route (1e-2 .. 1e-1) was
    off by three orders -- the tank shorts the baseband VOLTAGE, but the
    phase responds to the CURRENT through the PPV.  Anchored at the source
    (docs session, Hajimiri & Lee 1998): the ISF is defined per injected
    CHARGE ("amount of excess phase proportional to the ratio of the
    injected charge"; Fig. 6 is phase shift versus injected charge), so a
    per-charge sensitivity cannot be priced through a per-volt impedance --
    the units do not meet.  And the confirming shape here -- break the
    symmetry, check the lift is LINEAR in the breaking parameter -- is the
    paper's own linearity verification ("injecting impulses with different
    areas"), reproduced.
    """
    circuit.default_toolkit = circuit.numeric
    got = {}
    for a in (0.0, 0.05, 0.25):
        with quiet():
            cir, pss, pac = _lc_osc(a=a, psd=1e-6, npts=240)
            f0 = 1.0 / float(pss.period)
            _m_am, m_pm = pac.am_pm(pss, 1e-3 * f0, 0, harmonic=1)
        got[a] = float(np.abs(np.asarray(m_pm)).max())
    assert got[0.05] / got[0.0] > 1e6, got
    assert abs((got[0.25] / got[0.05]) / 5.0 - 1.0) < 0.02, got


def test_pac_solve_and_the_adjoint_transfer_row_take_the_deflated_route_on_an_oscillator():
    """2026-09-08: `_deflated_solve` (the Gourary removal, bordered with
    both null vectors) was wired into `adjoint_sideband_row` only;
    `PAC.solve` and `adjoint_transfer_row` solved plain outside
    `HARMONIC_GUARD`, where the plain solve's relative error is
    `eta / (2 pi df/f0)` with `eta = |lambda_1 - 1|` (docs session, four
    digits over five decades).  Now all three take the deflated route on an
    autonomous circuit and the plain one on a driven circuit.  Pinned: on
    the van der Pol under gear, with offsets r from the CARRIER (near DC
    the pole is not excited: the symmetric orbit's PPV has no DC term),
    (1) the deflated sweep equals the plain solve where the plain one is
    trustworthy (r = 1e-3, to 1e-6; measured 1.2e-8, the plain GMRES
    tolerance); (2) the deflated answer carries the physical pole, |y|
    scaling as 1/r between r = 1e-9 and 1e-10 to 1 %; (3) at r = 1e-10 the plain solve either
    refuses (measured: GMRES residual 1.7e-6 on the near-singular operator)
    or differs from it by more than 1e-6 relative; (4) `deflated` is True there and False on
    the driven mixer, whose result is unchanged.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0); cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    cir['ac'] = IS('v', gnd, i=0.0, iac=1.0)          # the sweep's small-signal source
    pss = PSS(cir, method='gear', reltol=1e-12)
    with quiet(AccuracyWarning):
        pss.solve(period=6.66, timestep=6.66 / 240, x0=np.array([2.0, 0.0]), maxiterations=80)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / float(pss.period)
    fp = pss.factored_period()
    T = float(fp.T)
    with quiet():
        ## the source vector for the sweep, as `solve` forms it (the AC
        ## vector with the reference row removed)
        u_ac = np.delete(np.asarray(cir.u(0.0, analysis=pac.par.analysis), dtype=complex).ravel(),
                         pss.irefnode)
        def rhs_and_alpha(r):
            f = r * f0
            w, _ = pss._forced_replay(fp, f, u_ac)
            a = np.exp(-2j * np.pi * f * T)
            return a, a * np.asarray(w)
        tol = max(pss.par.reltol * pac.KRYLOV_FACTOR, 1e-14)
        ## ⚠ NEAR THE CARRIER, not near DC: the pole at k = 0 is not
        ## excited on this half-wave-symmetric orbit (the PPV has no DC
        ## coefficient), so offsets are taken from f0, where Gamma_1 != 0
        a, b = rhs_and_alpha(1.0 + 1e-3)
        y_plain = pac._solve_each(fp, [a], [b], tol)[0][0]
        y_defl = pac._deflated_solve(pss, a, b, transposed=False, tol=tol)
        assert np.linalg.norm(y_defl - y_plain) < 1e-6 * np.linalg.norm(y_plain)   # measured 1.2e-8: the plain GMRES tolerance
        ## the sweep itself now takes that route
        res = pac.solve(pss, [(1.0 + 1e-3) * f0], sweeptype='absolute')
        assert res.info['deflated'] is True and res.info['matvecs'] is None
        ## the pole, carried analytically
        a9, b9 = rhs_and_alpha(1.0 + 1e-9); a10, b10 = rhs_and_alpha(1.0 + 1e-10)
        y9 = pac._deflated_solve(pss, a9, b9, transposed=False, tol=tol)
        y10 = pac._deflated_solve(pss, a10, b10, transposed=False, tol=tol)
        assert abs(np.linalg.norm(y10) / np.linalg.norm(y9) / 10.0 - 1.0) < 0.01
        ## the plain solve there: measured, it REFUSES (GMRES residual
        ## 1.7e-6 against a near-singular operator), which is the stronger
        ## form of "off by eta/(2 pi r)"; either outcome is the point
        try:
            y10_plain = pac._solve_each(fp, [a10], [b10], tol)[0][0]
        except RuntimeError as e:
            assert 'near-singular' in str(e) or 'did not converge' in str(e), str(e)
        else:
            assert np.linalg.norm(y10_plain - y10) > 1e-6 * np.linalg.norm(y10), \
                'the plain solve should be off by eta/(2 pi r) here'
        ## the adjoint row takes it too
        _row, row_info = pac.adjoint_transfer_row(pss, (1.0 + 1e-3) * f0, 0)
        assert row_info['deflated'] is True
    ## driven: plain route, unchanged
    from pycircuit.circuit.elements import Diode
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    c['R'] = R(1, 2, r=1e4); c['D'] = Diode(2, gnd); c['C'] = C(2, gnd, c=1e-12)
    p2 = PSS(c, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        p2.solve(period=1e-6, timestep=1e-6 / 200, maxiterations=40)
        pac2 = PAC(c, toolkit=circuit.numeric)
        res2 = pac2.solve(p2, [0.13e6])
        assert res2.info['deflated'] is False and res2.info['matvecs'] is not None
        _row, row_info = pac2.adjoint_transfer_row(p2, 0.13e6, 2)
        assert row_info['deflated'] is False


def test_the_forward_pac_sidebands_equal_the_adjoint_rows_on_a_non_uniform_grid():
    """The forward `PAC.solve` and `adjoint_sideband_row` must use the SAME
    period weights -- on a 3:1 grid too.  The uniform-grid twin of this check
    is `test_pac_reports_sidebands_at_the_right_frequencies_and_conjugates_
    the_fold`; this one runs the trapezoid-weighted branch of both sides.
    ⚠ It pins CONSISTENCY only: the identity closes for any weights used on
    both sides, which is why it never caught the old `1/N`.  Correctness on a
    non-uniform grid is pinned by
    `test_pnoise_is_correct_on_a_non_uniform_grid_with_trapezoid_period_weights`.
    """
    circuit.default_toolkit = circuit.numeric
    fclk = 100e3; T = 1.0 / fclk; fin = 7e3; n = 200
    cir = SubCircuit()
    for nd in ('in', 'out', 'ck'):
        cir.add_node(nd)
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0, vac=1.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9,
                           vth=0.0, vs=50e-3)
    cir['C0'] = C('out', gnd, c=100e-12)
    w = 1.0 + 0.5 * np.sin(2.0 * np.pi * np.arange(n) / n)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / n, x0=np.zeros(cir.n - 1),
                  grid=w / w.sum(), maxiterations=60, break_events=False)
    assert pss.converged
    assert pss._period_quadrature(pss.factored_period()) is not None
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(pss, freqs=[fin])
    fout = np.asarray(res.sweep_values, dtype=float)
    X = np.asarray(res.x)
    io = [str(x) for x in cir.nodes].index('out')
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), pss.irefnode,
                             circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    H = np.asarray(pac.adjoint_sideband_row(pss, fin, io, sidebands=[0, 1, -1])[0])
    for li, l in enumerate((0, 1, -1)):
        f_phys = abs(fin + l * fclk)
        k = int(np.argmin(np.abs(fout - f_phys)))
        assert abs(fout[k] - f_phys) < 1e-6 * fclk, (l, fout[k], f_phys)
        x = complex(X[io, k]); h = complex(H[li] @ u_ac)
        if fin + l * fclk < 0:
            h = np.conj(h)
        assert abs(x - h) < 1e-10 * abs(h), (l, x, h)


def test_one_step_factored_periods_replay_on_the_grid_the_solve_was_on():
    """The factored period of a stage method (radau, trbdf2, esdirk43, and
    the TR-BDF2 twin trap/euler borrow) replays the converged orbit on the
    grid the SOLVE ran on, not on a uniform `linspace` of the same count
    (2026-09-20, Andreas: "honour the caller's grid").

    ⚠ UNTIL THEN EVERY ONE-STEP REPLAY WAS UNIFORM, and the consequence
    was live on the DEFAULT method: `break_events` lands a pulse edge in
    the solve, `factored_period_full` replayed on a uniform grid, so
    `_period_quadrature(fp)` returned None and `carrier_phasor` took a
    plain MEAN over the event grid's non-uniform nodes.  Pulsed RC, tr = 0,
    radau, events landed by default, first harmonic of the capacitor
    voltage against the analytic periodic solution's coefficient::

        N      nodes   uniform replay (shipped)   replay on the solved grid
        50     55      1.04e-01                   1.54e-03
        100    105     5.43e-02                   3.77e-04
        200    205     2.77e-02                   9.39e-05
        400    403     7.05e-03                   2.33e-05
        800    803     3.53e-03                   5.68e-06

    First order against second, 600x at 800.  (2026-09-21: with the period
    quadrature breaking its spline at the landed nodes the same three
    grids read 7.9e-8 / 4.6e-8 / 1.4e-9 -- the reference's floor.)  The
    same uniform replay is
    why radau's PPV, modes and noise read as "exact on the 3:1 grid": its
    adjoint never saw that grid.  `factored_period` now hands the solved
    fractions down (`_replay_grid`); a direct `factored_period_stage(x0, T,
    npts)` call with a bare count is uniform as before, and a uniform
    solve is bit-identical (the fractions are only passed when the grid is
    not uniform).
    """
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    RR, CC = 1e3, 3e-10
    TAU = RR * CC
    TD, PW = 0.0125 * T, 0.4 * T

    def pulsed():
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=TD, tr=0.0, tf=0.0, pw=PW, per=T)
        c['R'] = R(1, 2, r=RR)
        c['C'] = C(2, gnd, c=CC)
        return c

    def exact(ts):
        tr = 1e-18
        e = [0.0, TD, TD + tr, TD + tr + PW, TD + tr + PW + tr, T]
        seg = [(e[0], e[1], 0.0, 0.0), (e[1], e[2], 0.0, 1.0 / tr),
               (e[2], e[3], 1.0, 0.0), (e[3], e[4], 1.0, -1.0 / tr),
               (e[4], e[5], 0.0, 0.0)]

        def prop(x0, t0, t1, a, b):
            one_e = -np.expm1(-(t1 - t0) / TAU)
            return x0 * (1.0 - one_e) + (a - b * TAU) * one_e + b * (t1 - t0)
        alpha, beta = 1.0, 0.0
        for (t0, t1, a, b) in seg:
            al = np.exp(-(t1 - t0) / TAU)
            alpha, beta = al * alpha, al * beta + prop(0.0, t0, t1, a, b)
        x0 = beta / (1.0 - alpha)
        out = np.empty(len(ts))
        for i, t in enumerate(ts):
            t = float(t) % T
            x = x0
            for (t0, t1, a, b) in seg:
                if t <= t1 + 1e-30:
                    out[i] = prop(x, t0, t, a, b)
                    break
                x = prop(x, t0, t1, a, b)
            else:
                out[i] = x
        return out
    tt = np.linspace(0.0, T, 100001)[:-1]
    X1 = complex(np.sum(exact(tt) * np.exp(-2j * np.pi * tt / T)) / len(tt))

    errs = []
    for N in (100, 200, 400):
        c = pulsed()
        p = PSS(c, method='radau', reltol=1e-10)
        with quiet():
            p.solve(period=T, timestep=T / N, maxiterations=40)
        assert p.converged and p.break_events and len(p.event_times) == 4
        fp = p.factored_period()
        ## the replay IS on the solved grid
        assert np.allclose(np.asarray(fp.times, float),
                           np.asarray(p.waveform[0], float).ravel(), rtol=0, atol=1e-18 * 0 + 1e-15)
        assert p._period_quadrature(fp) is not None
        full = [str(n_) for n_ in c.nodes].index('2')
        out = full if full < p.irefnode else full - 1
        a = PAC(c).carrier_phasor(p, out)
        errs.append(abs(a - X1) / abs(X1))
    assert errs[0] < 1e-3 and errs[2] < 5e-5, errs         # shipped: 5.4e-2 / 7.0e-3
    ## and on this EVENT grid the weights are the PIECEWISE spline's,
    ## breaking at the landed nodes (2026-09-21, item 2 of the
    ## non-uniform-grid list; until then the trapezoid, "a spline rings
    ## through the kink" -- true of a spline that crosses it).  Measured
    ## here against the analytic coefficient: 3.77e-4 / 9.39e-5 / 2.33e-5
    ## under the trapezoid (second order), 7.9e-8 / 4.6e-8 / 1.4e-9 under
    ## the piecewise rule -- at the reference's own floor, hence no ladder
    _tm = np.asarray(fp.times, float)
    _N = len(fp.steps)
    _h = np.diff(_tm[:_N + 1])
    assert p.event_times
    from pycircuit.circuit.shooting import periodic_spline_weights as _psw
    _Tq = _tm[_N] - _tm[0]
    np.testing.assert_allclose(p._period_quadrature(fp) * _Tq,
                               _psw(_tm[:_N], _Tq, p._event_nodes(_tm[:_N], _Tq)),
                               rtol=1e-12, atol=0)
    assert np.max(np.abs(p._period_quadrature(fp) * _Tq - 0.5 * (_h + np.roll(_h, 1)))) > 0.1 * np.max(_h)
    assert errs[0] < 5e-7 and errs[2] < 2e-8, errs

    ## a bare count is still a uniform replay; a uniform solve is unchanged
    c = pulsed()
    p = PSS(c, method='radau', reltol=1e-10)
    with quiet():
        p.solve(period=T, timestep=T / 100, maxiterations=40, break_events=False)
    fpu = p.factored_period()
    assert p._period_quadrature(fpu) is None
    x0 = np.asarray(p._period_state[1], float).ravel()
    fpb = p.factored_period_stage(x0, T, 100)
    assert np.allclose(np.asarray(fpb.times, float), np.linspace(0.0, T, 101), rtol=0, atol=1e-22)

@pytest.mark.parametrize('method,kind', [('radau', 'stage'), ('glm2', 'glm'),
                                         ('glm3', 'glm')])
def test_pac_on_a_staged_solve_borders_its_sideband_solve_with_the_event_rows_and_is_exact(method, kind):
    """Phase B of events-as-unknowns (2026-09-22): on a solve whose grid was
    landed on state events, a periodic perturbation moves the crossings,
    and `PAC.solve` borders its per-frequency system with the event rows
    (`_event_columns`; block elimination, a K x K Schur complement for the
    crossings' modulation).  Verified against the finite difference AT
    FIXED BASE -- an inner Newton on `(x_0, theta)` with the perturbed
    circuits on the solve's own grid, the discrete map the bordering
    linearises: PWM loop, 100 points, `out` 4e-10, `fb` 1.4e-10, `sw`
    4e-8 at f0 and 2 f0, the crossings' shifts identical to six digits;
    unbordered 60-70 % off on the filter nodes and 3-12x on the switch
    node.  Read from `PAC.time_response`: the sideband coefficients are
    quadrature integrals of the envelope and on this grid (eight steps of
    6e-4 T inside the switch's window) summing them back at the nodes is
    not exact.  ⚠ Two instrument traps, recorded because they cost a day: `vac`
    DEFAULTS TO 1 on every voltage source, so every source in a fixture is
    an AC source unless zeroed; and re-SOLVED orbits at +-eps each re-land
    their base grid from their own first stage, so their difference is
    the derivative of a different discrete map (9 % off here) -- the FD
    must keep the base.  Pinned at 60 points, f = f0: the reconstruction
    of `PAC.result` (AC-phasor convention, divided by the source's AC
    phase) against the fixed-base FD to 1e-6 on `out` and `fb`, the
    unbordered one at least 0.1 / 0.03 off (0.55 / 0.057 measured), and
    `event_shifts` against the FD.  2026-09-25: the Nordsieck GLMs on their
    own maps (the fixed-base walk is theirs, `_walk('glm')`): 6e-10 / 7e-10
    (glm2), 8.7e-10 / 6.5e-10 (glm3), the shifts 7e-10 / 6e-10, the
    unbordered response 3e-3 / 5e-3 .. 7e-3 off.
    """
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    f0 = 1.0 / T
    N = 60

    def build(eps=0.0, vac=0.0):
        cir = _pwm_loop(T)
        del cir['Vin']
        cir.add_node('vin0')
        cir['Vin'] = VS('vin0', gnd, v=5.0, vac=0.0)
        cir['Vp'] = VSin('vin', 'vin0', vo=0.0, va=eps, freq=f0, phase=90.0, vac=vac)
        cir['Vramp'].iparv.vac = 0.0
        return cir

    cir = build(0.0, vac=1.0)
    p0 = PSS(cir, method=method, reltol=1e-10)
    with quiet():
        p0.solve(period=T, timestep=T / N, x0=np.zeros(cir.n - 1), maxiterations=100)
    assert p0.converged and p0._event_columns is not None
    g0 = np.asarray(p0._grid_fracs, float)
    x0s = np.asarray(p0._period_state[1], float)
    m = cir.n - 1
    th0 = np.asarray(p0._state_event_fracs, float)
    K = len(th0)
    ev = p0._event_columns
    Wk, ck = np.asarray(ev['W']), np.asarray(ev['c'])

    def solve_fixed_base(eps):
        c2 = build(eps)
        q = PSS(c2, method=method, reltol=1e-10)
        with quiet(ConvergenceWarning):
            q.solve(period=T, timestep=T / len(g0), x0=x0s, grid=g0,
                    maxiterations=1, state_events=False)
        z = np.concatenate((x0s, th0))
        for _it in range(40):
            xx, th = z[:m], z[m:]
            fr, hsens, nodes = q._event_remap(g0, th0, th, T)
            hs = fr * T
            tms = np.concatenate(([0.0], np.cumsum(hs)))
            wk = q._walk(kind, xx, tms, hs, T=T, hsens=hsens,
                         capture=set(range(1, len(hs) + 1)))
            x_end, Mx, Pk = wk.x_end, wk.P, wk.Pk
            F = np.zeros(m + K)
            J = np.zeros((m + K, m + K))
            F[:m] = xx - np.asarray(x_end)
            J[:m, :m] = np.eye(m) - np.asarray(Mx)
            for k in range(K):
                J[:m, m + k] = -np.asarray(Pk[k]).ravel()
            for k, nd in enumerate(nodes):
                xj, Pj, Pkj = q._captured[nd]
                F[m + k] = float(Wk[k] @ xj) - ck[k]
                J[m + k, :m] = Wk[k] @ Pj
                for l in range(K):
                    J[m + k, m + l] = float(Wk[k] @ np.asarray(Pkj[l]).ravel())
            if np.max(np.abs(F)) < 1e-12:
                break
            z = z - np.linalg.solve(J, F)
        X = np.array([z[:m]] + [np.asarray(q._captured[j][0]) for j in range(1, len(hs) + 1)]).T
        return z, X, tms

    e = 1e-4
    zp, Xp, tms = solve_fixed_base(e)
    zm, Xm, tms_m = solve_fixed_base(-e)
    d = (Xp - Xm) / (2 * e)
    dth_fd = (zp[m:] - zm[m:]) / (2 * e)
    ## the FD reads "node j" of the two solves, whose TIME moved with the
    ## crossings; the response PAC reports is at FIXED times (2026-09-22),
    ## so take the node's motion out with the orbit's rate -- the
    ## fixed-time form itself is exact-verified on the autonomous fixture
    ## (`test_the_sideband_response_on_a_staged_oscillator_is_bordered...`)
    _rate = PAC(cir, toolkit=circuit.numeric)._orbit_rate(p0, p0._event_columns['nodes'])   # (N+1, m)
    d = d - _rate.T[:, :d.shape[1]] * ((tms - tms_m)[:d.shape[1]] / (2 * e))[None, :]
    names = [str(n_) for i, n_ in enumerate(cir.nodes) if i != p0.irefnode]
    errs = {}
    for bordered in (True, False):
        if not bordered:
            p0._event_columns = None
        pac = PAC(cir, toolkit=circuit.numeric)
        with quiet():
            res = pac.solve(p0, [f0])
        t_r, y_r = res.info['time_response'][0]
        assert np.max(np.abs(t_r - tms[:len(t_r)])) < 1e-5 * T      # the same nodes (fp rebuilds them from the fractions)
        rec = np.real(y_r / 1j).T                                          # the source's AC phase is 90 deg
        errs[bordered] = {nm: np.max(np.abs(rec[names.index(nm)] - d[names.index(nm)][:rec.shape[1]]))
                          / np.max(np.abs(d[names.index(nm)])) for nm in ('out', 'fb')}
        if bordered:
            sh = np.real(np.asarray(res.info['event_shifts'][0]) / 1j)
            assert np.max(np.abs(sh - dth_fd)) < 1e-6 * np.max(np.abs(dth_fd)), (sh, dth_fd)
    p0._event_columns = ev
    assert errs[True]['out'] < 1e-6 and errs[True]['fb'] < 1e-6, errs
    ## ⚠ RE-PINNED 2026-09-22 with the fixed-time instrument: the unbordered
    ## response is 3.1 % / 4.9 % off at 60 points (out / fb), NOT the 55 % /
    ## 5.7 % first recorded -- that 55 % was mostly the moving-node artefact
    ## the old instrument shared with the old `time_response`.  The
    ## bordering's gain on this driven loop is real but modest; the O(1)
    ## cases are the staged oscillator's sideband response and the
    ## covariance closure (their own tests).
    ## ⚠ RE-PINNED AGAIN 2026-09-22 for VSwitch's COMPACT transition: 1.1 % /
    ## 2.5 % (was 3.1 % / 4.9 % with the tanh, whose tails outside the window
    ## a fixed grid never resolved).  The bordering is the exact
    ## linearisation; its numerical weight on this loop is a few per cent.
    ## ⚠ RE-PINNED A THIRD TIME 2026-09-25, for the 16-step window
    ## (`event_window_steps`): 0.076 % / 0.16 % -- the rest of the old gap
    ## was the 8-step window's own resolution
    assert errs[False]['out'] > 5e-4 and errs[False]['fb'] > 1e-3, errs     # 7.6e-4 / 1.6e-3 at 60 points


@pytest.mark.parametrize('method', ['radau', 'glm2', 'glm3', 'trap', 'euler'])
def test_the_adjoint_sideband_row_on_a_staged_solve_is_the_transpose_of_the_bordered_forward_solve(method):
    """Phase B of events-as-unknowns (2026-09-22): `adjoint_sideband_row`
    borders its adjoint solve with the event rows -- the transpose of
    `PAC.solve`'s bordered system, by block elimination and a second
    reverse pass carrying the event rows' costate `-zeta_k w_k` -- so
    `pnoise` and `mixer_response` on a solve whose grid was landed on
    state events read the moving crossings.  Pinned as the suite pins
    PAC everywhere: dual consistency, the bordered forward solve's
    reported coefficients equal `H_l . u_ac` to 1e-10 for sidebands
    0, 1, -1 (PWM loop, 60 points, f = 0.3 f0), and the unbordered row
    differs from them by more than 1 %.  2026-09-25: the Nordsieck GLMs on
    their own maps too (two growth restarts in this period): 1e-15, the
    unbordered row 1.1e-3 .. 1.9e-3 off (radau 3.5e-4 .. 4.6e-4).
    2026-09-30 (the review, F7): the PLAIN map too -- trap and euler were
    never admitted to the bordering and missed the forward solve by
    1.9e-3 / 9e-3; now 3.6e-13 / 4e-15.
    """
    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    f0 = 1.0 / T
    fin = 0.3 * f0
    cir = _pwm_loop(T)
    del cir['Vin']
    cir.add_node('vin0')
    cir['Vin'] = VS('vin0', gnd, v=5.0, vac=0.0)
    cir['Vp'] = VSin('vin', 'vin0', vo=0.0, va=0.0, freq=fin, phase=0.0, vac=1.0)
    cir['Vramp'].iparv.vac = 0.0
    p = PSS(cir, method=method, reltol=1e-10)
    with quiet(AccuracyWarning, ConvergenceWarning):
        p.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1), maxiterations=100)
    assert p.converged and p._event_columns is not None
    io_full = [str(n_) for n_ in cir.nodes].index('out')
    io = io_full if io_full < p.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), p.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(p, freqs=[fin])
    fout = np.asarray(res.sweep_values, dtype=float)
    X = np.asarray(res.x)
    H = np.asarray(pac.adjoint_sideband_row(p, fin, io, sidebands=[0, 1, -1])[0])
    ev = p._event_columns
    p._event_columns = None
    H0 = np.asarray(pac.adjoint_sideband_row(p, fin, io, sidebands=[0, 1, -1])[0])
    p._event_columns = ev
    for li, l in enumerate((0, 1, -1)):
        f_phys = abs(fin + l * f0)
        k = int(np.argmin(np.abs(fout - f_phys)))
        x = complex(X[io_full, k])
        h = complex(H[li] @ u_ac)
        h0 = complex(H0[li] @ u_ac)
        if fin + l * f0 < 0:
            h = np.conj(h)
            h0 = np.conj(h0)
        assert abs(x - h) < 1e-10 * abs(h), (l, x, h)
        ## the unbordered row: 0.55 % off with VSwitch's compact transition
        ## (2026-09-22; several per cent with the tanh's tails); 0.036 % with
        ## the window in 16 steps (2026-09-25)
        assert abs(x - h0) > 1e-4 * abs(h), (l, x, h0)

@pytest.mark.parametrize('method', ['radau', 'gear', 'trap'])
def test_the_adjoint_transfer_row_on_a_staged_solve_is_the_transpose_of_the_bordered_forward_solve(method):
    """The review's D3 (2026-10-01): `adjoint_transfer_row` -- every source
    to the output at `t = 0` -- was never bordered on a staged solve, while
    `PAC.solve` and `adjoint_sideband_row` are: on the PWM loop (60 points,
    f = 0.3 f0) it missed the bordered forward solve's state at `t = 0` by
    2.8e-4 (radau) / 1.6e-3 (trap) / 3.9e-2 (gear), exactly the unbordered
    row's gap.  Bordered as the sideband family is (the event rows'
    costate, a second reverse pass): 2.7e-15 / 1.0e-13 / 6.8e-15."""

    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    fin = 0.3 / T
    cir = _pwm_loop(T)
    del cir['Vin']
    cir.add_node('vin0')
    cir['Vin'] = VS('vin0', gnd, v=5.0, vac=0.0)
    cir['Vp'] = VSin('vin', 'vin0', vo=0.0, va=0.0, freq=fin, phase=0.0, vac=1.0)
    cir['Vramp'].iparv.vac = 0.0
    p = PSS(cir, method=method, reltol=1e-10)
    with quiet(AccuracyWarning, ConvergenceWarning):
        p.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1), maxiterations=100)
    assert p.converged and p._event_columns is not None
    io_full = [str(n_) for n_ in cir.nodes].index('out')
    io = io_full if io_full < p.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), p.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(p, freqs=[fin], sweeptype='absolute')
        y0 = complex(np.asarray(res.info['time_response'][0][1])[0][io])
        h = complex(pac.adjoint_transfer_row(p, fin, io)[0] @ u_ac)
        ev = p._event_columns
        p._event_columns = None
        h0 = complex(pac.adjoint_transfer_row(p, fin, io)[0] @ u_ac)
        p._event_columns = ev
    assert abs(h - y0) < 1e-11 * abs(y0), (h, y0)
    assert abs(h0 - y0) > 1e-4 * abs(y0), (h0, y0)


def test_the_oscillator_consumers_run_bordered_on_a_staged_gear_solve():
    """Gear's staged OSCILLATOR (the autonomous pair stage, 2026-09-24) was a
    combination no consumer had met, and three broke on it:

    * `floquet_modes` CRASHED: the event costate injection was built at the
      MAP's width (`2m` on the pair) and the pair's reverse step adds it to
      the `m`-wide circuit block;
    * the bordered `PAC.solve` CRASHED twice over: the source's forced
      replay seeded at `m`, and `dtheta/dz` (K x 2m) contracted with the
      first `m` of the response;
    * the ADJOINT row (pnoise's) stayed UNBORDERED on an autonomous pair --
      a deliberate exclusion while no such solve existed -- and missed the
      bordered forward solve by 8.4 % / 2.1 % (sidebands 0 / 1).

    Measured after: the forced response 1.2e-3 from the exact one at 200
    points (5.2e-5 at 800) against 2.2e-2 unbordered, the adjoint row the
    transpose of the forward solve to 5e-15, the PPV 6.3e-3 from the exact
    saltation PPV, and `diffusion_constant` / `pnoise` / `oscillator_
    covariance` within 0.04 % / 3 % / 0.6 % of radau's at 200 points (0.00 /
    0.1 / 0.3 % at 800) -- gear's own second order.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    cir['Vdd'] = VS('vdd', gnd, v=5.0, vac=1.0)
    names = [str(n_) for n_ in cir.nodes]
    seed, Tl = _relaxation_oscillator_seed(cir)
    q = PSS(cir, method='gear', reltol=1e-9)
    with quiet(AccuracyWarning):
        q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100,
                state_events=True)
    assert q.converged and q._event_columns is not None
    Tq = float(q.period)
    sel = [names.index(nm) for nm in ('c', 'fb0', 'fb1')]
    red = [n_ for i, n_ in enumerate(names) if i != q.irefnode]
    idx = [red.index(nm) for nm in ('c', 'fb0', 'fb1')]
    X = np.asarray(q.waveform[1], dtype=float)
    ts = np.asarray(q.waveform[0], dtype=float)

    ## the Floquet modes run, on the TOTAL map
    fp = q.factored_period()
    Md = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                          for e in np.eye(fp.width)])
    Mt = Md + (np.asarray(q._event_columns['P_end'], dtype=float)
               @ np.asarray(q._event_sensitivity, dtype=float))
    lam_t = np.sort(np.abs(np.linalg.eigvals(Mt)))[::-1]
    with quiet(UsageWarning):
        fm = q.floquet_modes(nmodes=2)
        v, info = q.ppv()
    lam_fm = np.sort(np.abs([mm['lam'] for mm in fm]))[::-1]
    assert np.max(np.abs(lam_fm - lam_t[:2]) / lam_t[:2]) < 1e-9, (lam_fm, lam_t[:2])

    ## the PPV against the exact saltation PPV, before and after the crossings
    _T_ex, exact = _exact_relaxation_oscillator_ppv()
    S = np.asarray(info['samples'])
    nodes = list(q._event_columns['nodes'])
    js = ([j for j in (5, 15, 25) if j < nodes[0]]
          + [int(np.searchsorted(ts, f * Tq)) for f in (0.5, 0.8)])
    worst = max(float(np.max(np.abs(S[j, idx] - exact(X[sel, j])))
                      / np.max(np.abs(exact(X[sel, j])))) for j in js)
    assert worst < 2e-2, worst

    ## the forced response against the exact one: bordered, and not
    T_ex2, _a0, orbit_at = _exact_relaxation_oscillator_forced(0.3, 0.0)
    tgrid = np.linspace(0.0, T_ex2, 20001)[:-1]
    orb = np.array([orbit_at(tt) for tt in tgrid])
    t_shift = float(tgrid[int(np.argmin(np.linalg.norm(orb - X[sel, 0], axis=1)))])
    _T, at, _o = _exact_relaxation_oscillator_forced(0.3, t_shift)
    scale = np.max([np.abs(at((t_shift + float(x)) % T_ex2))
                    for x in np.linspace(0.0, Tq, 400)], axis=0)
    probe = [0, 20, int(np.searchsorted(ts, 0.35 * Tq)), int(np.searchsorted(ts, 0.7 * Tq))]
    ev = q._event_columns
    worst = {}
    for bordered in (True, False):
        q._event_columns = ev if bordered else None
        try:
            pac = PAC(cir, toolkit=circuit.numeric)
            with quiet():
                res = pac.solve(q, [0.3 / Tq], sweeptype='absolute')
        finally:
            q._event_columns = ev
        tt, yy = res.info['time_response'][0]
        worst[bordered] = max(
            float(np.max(np.abs(np.asarray(yy[j])[idx]
                                - at((t_shift + float(tt[j])) % T_ex2)
                                * np.exp(1j * 2.0 * np.pi * 0.3
                                         * ((t_shift + float(tt[j])) // T_ex2)))
                         / scale)) for j in probe)
    ## ⚠ RE-PINNED 2026-09-25 for the 16-step window: bordered 1.1e-3,
    ## unbordered 4.8e-3 (it was above 1e-2 at 8 steps)
    assert worst[True] < 5e-3 and worst[False] > 3e-3, worst

    ## the adjoint row is the transpose of the same bordered solve
    fin = 0.3 / Tq
    io_full = names.index('fb1')
    io = io_full if io_full < q.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), q.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(q, freqs=[fin], sweeptype='absolute')
        H = np.asarray(pac.adjoint_sideband_row(q, fin, io, sidebands=[0, 1])[0])
    fout = np.asarray(res.sweep_values, dtype=float)
    Xr = np.asarray(res.x)
    for li, l in enumerate((0, 1)):
        k = int(np.argmin(np.abs(fout - (fin + l / Tq))))
        x = complex(Xr[io_full, k])
        h = complex(H[li] @ u_ac)
        assert abs(x - h) < 1e-10 * abs(h), (l, x, h)


def _exact_relaxation_oscillator_forced(fr, t_shift):
    """The EXACT sideband response of `_comparator_relaxation_oscillator`
    to a unit tone on the rail at `fr` times its own fundamental, on
    `_exact_relaxation_oscillator_model`: the variational system ``y' =
    A_s y + B e^{j w (t - t_shift)}`` in each switch state (augmented
    matrix exponentials), ``y+ = S y-`` at the crossings, and ``y(T) =
    e^{j w T} y(0)``.  `t_shift` is where the PSS's `t = 0` sits on this
    orbit.  Returns ``(T, at, orbit_at)`` with ``at(t)`` the response and
    ``orbit_at(t)`` the state, `t` from the ON->OFF switching."""
    from scipy.linalg import expm
    mdl = _exact_relaxation_oscillator_model()
    w = 2.0 * np.pi * fr / mdl.T

    def aug(A, t):
        Z = np.zeros((4, 4), dtype=complex)
        Z[:3, :3] = A
        Z[:3, 3] = mdl.Bv
        Z[3, 3] = 1j * w
        return expm(Z * t)

    def seg(A, y_in, ph_in, t):
        Z = aug(A, t) @ np.concatenate((y_in, [ph_in]))
        return Z[:3], Z[3]

    ph0 = np.exp(-1j * w * t_shift)
    y, ph = seg(mdl.A_off, mdl.S0 @ np.zeros(3, dtype=complex), ph0, mdl.t_off)
    y, ph = seg(mdl.A_on, mdl.S1 @ y, ph, mdl.t_on)
    alpha = np.exp(-1j * w * mdl.T)
    y0 = np.linalg.solve(np.eye(3) - alpha * mdl.M, alpha * y)

    def at(t):
        yy, pp = mdl.S0 @ y0, ph0
        if t < mdl.t_off:
            yy, pp = seg(mdl.A_off, yy, pp, t)
        else:
            yy, pp = seg(mdl.A_off, yy, pp, mdl.t_off)
            yy, pp = seg(mdl.A_on, mdl.S1 @ yy, pp, t - mdl.t_off)
        return yy

    return mdl.T, at, mdl.orbit_at


def test_the_sideband_response_on_a_staged_oscillator_is_bordered_deflated_and_matches_the_exact_forced_response():
    """Events phase B, last item (2026-09-22): PAC on an AUTONOMOUS staged
    solve.  The bordered system collapses onto the total map -- `(I - a
    M_tot) y_0 = a (w + P_theta dtheta_f)`, `dtheta_f = -Gt^-1 W f_node`
    the source's own motion of the crossings -- solved by the deflated
    route with the TOTAL operator (the fixed-grid `M` has no unit
    multiplier), and the nodes' responses are at FIXED times (`Pk - xdot
    tau^T`).  The reference is EXACT (`_exact_relaxation_oscillator_forced`)
    and the comparison sits at the same offset from each model's own f0.
    Measured at 200 points with VSwitch's compact transition
    (2026-09-23, each component scaled by its own orbit maximum): 9.8e-5 /
    2.5e-4 / 1.5e-4 at 0.3 / 1.7 / 1.001 f0, every node and component
    included -- the collapsed `c` inside the 10 ns ON phase since the rate
    became the DAE's own (the three-node stencil had cost it 5 %), and
    `fb0` at its near-null at t/T = 0.35, which a per-node relative error
    had been reading as 7.4e-3 while its absolute error was the same 1e-5
    as everywhere on the orbit.  The PLAIN deflated
    solve on the same staged solve reads the same to 4 digits: with the whole
    transition inside the landed window the fixed-grid map carries the
    switching.  With the tanh it read 0.3-400x off and the bordered one
    15.7 % at 1.001 f0 (E8): the tails outside the window.  The adjoint row is the transpose of the same solve:
    dual-consistent with the forward one (the deflated solve refines on
    the plain operator, so both are the discrete operator's own)."""
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    cir['Vdd'] = VS('vdd', gnd, v=5.0, vac=1.0)
    names = [str(n_) for n_ in cir.nodes]
    seed, Tl = _relaxation_oscillator_seed(cir)
    q = PSS(cir, method='radau', reltol=1e-9)
    with quiet():
        q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=100, state_events=True)
    assert q.converged and q._event_columns is not None
    Tq = float(q.period)
    ts = np.asarray(q.waveform[0], dtype=float)
    X = np.asarray(q.waveform[1], dtype=float)
    sel = [names.index(nm) for nm in ('c', 'fb0', 'fb1')]
    red = [n_ for i, n_ in enumerate(names) if i != q.irefnode]
    idx = [red.index(nm) for nm in ('c', 'fb0', 'fb1')]
    T_ex, _at0, orbit_at = _exact_relaxation_oscillator_forced(0.3, 0.0)
    tgrid = np.linspace(0.0, T_ex, 20001)[:-1]
    orb = np.array([orbit_at(t) for t in tgrid])
    k0 = int(np.argmin(np.linalg.norm(orb - X[sel, 0], axis=1)))
    t_shift = float(tgrid[k0])
    nodes = list(q._event_columns['nodes'])
    on_phase = set(range(nodes[1], nodes[2] + 1))          # between the two crossings: the switch is ON
    probe = [j for j in (0, 20, int(np.searchsorted(ts, 0.35 * Tq)), int(np.searchsorted(ts, 0.7 * Tq)))]
    ev = q._event_columns
    for fr in (0.3, 1.7, 1.001):
        _T, at, _o = _exact_relaxation_oscillator_forced(fr, t_shift)
        f = fr / Tq
        ## ⚠ EACH COMPONENT IS SCALED BY ITS OWN ORBIT MAXIMUM (2026-09-23),
        ## not by its value at the node.  `fb0`'s response passes through a
        ## near-null at t/T = 0.35 (1.3e-3 against its own 0.24), and a
        ## per-node relative error there measures the zero crossing, not the
        ## solve: that ONE entry read 7.4e-3 while its absolute error,
        ## 9.8e-6, was the same as everywhere else on the orbit, and it moved
        ## by 10x under changes that moved nothing else.  The trap is in the
        ## campaign's own notes ("a relative error against a zero-crossing
        ## quantity is not agreement") and this test had walked into it.
        scale = np.max([np.abs(at((t_shift + float(x)) % T_ex))
                        for x in np.linspace(0.0, Tq, 400)], axis=0)
        worst = {True: 0.0, False: 0.0}
        for bordered in (True, False):
            q._event_columns = ev if bordered else None
            try:
                pac = PAC(cir, toolkit=circuit.numeric)
                with quiet():
                    res = pac.solve(q, [f], sweeptype='absolute')
                tt, yy = res.info['time_response'][0]
            except RuntimeError:
                ## the plain deflated solve borders the fixed-grid map with
                ## the TOTAL map's null vectors: near a harmonic it does not
                ## even converge -- that is the defect showing, count it
                assert not bordered, (fr, 'the bordered deflated solve failed')
                worst[False] = np.inf
                continue
            finally:
                q._event_columns = ev
            for j in probe:
                te = (t_shift + float(tt[j])) % T_ex
                wrap = np.exp(1j * 2.0 * np.pi * fr * ((t_shift + float(tt[j])) // T_ex))
                ex = at(te) * wrap
                got = np.asarray(yy[j])[idx]
                ## every component at every probe, the collapsed `c` in the ON
                ## phase included: the 5 % it cost was the rate STENCIL, and the
                ## DAE derivative (2026-09-22) reads the exact rate to 1e-15
                worst[bordered] = max(worst[bordered], float(np.max(np.abs(got - ex) / scale)))
        ## ⚠ RE-PINNED 2026-09-22 for VSwitch's COMPACT transition: bordered
        ## 2.3e-3 / 7.4e-4 / 3.2e-3 and UNBORDERED 8.7e-4 / 6.5e-4 / 1.3e-3
        ## at 0.3 / 1.7 / 1.001 f0 -- both at the comparison's floor, the
        ## near-harmonic 15.7 % (E8: the tanh's multiplier displacement) gone.
        ## With the tanh the plain deflated solve read 0.3-400x off: the tails.
        ## ⚠ RE-PINNED 2026-09-23 on the component-scaled instrument:
        ## bordered 9.8e-5 / 2.5e-4 / 1.5e-4, unbordered the same to 4
        ## digits -- with the compact transition the landed grid's own map
        ## carries the switching, so the bordering is a sliver HERE (it is
        ## 10 % on gear's two-step map, and it is what makes the forward
        ## and adjoint solves one object).
        assert worst[True] < 1e-3, (fr, worst)
        assert worst[False] < 1e-3, (fr, worst)
    ## the adjoint row is the transpose of the same bordered, deflated solve
    fin = 0.3 / Tq
    io_full = names.index('fb1')
    io = io_full if io_full < q.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), q.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(q, freqs=[fin], sweeptype='absolute')
        H = np.asarray(pac.adjoint_sideband_row(q, fin, io, sidebands=[0, 1])[0])
    fout = np.asarray(res.sweep_values, dtype=float)
    Xr = np.asarray(res.x)
    for li, l in enumerate((0, 1)):
        k = int(np.argmin(np.abs(fout - (fin + l / Tq))))
        x = complex(Xr[io_full, k])
        h = complex(H[li] @ u_ac)
        assert abs(x - h) < 1e-8 * abs(h), (l, x, h)


def test_trap_opened_at_x0_is_the_transpose_of_its_forward_replay():
    """A trapezoidal plain map opened at `x(0)` (`x0_unknown=True`) starts
    with an order-dropped Euler step (``b = 0``) whose companion current
    `Pq` the next trap step (``b != 0``) reads.  `_LMMStep.adjoint` dropped
    the costate on `Pq` whenever ``b = 0`` -- right for gear, whose steps
    never read it -- so this map was not its forward replay's transpose.
    Found 2026-09-25 by the matrix-free event stage (its reverse-replayed
    event rows read 8.4e-4 off the dense ones).  Measured before the fix:
    the map's own duality 2.8e-3 on the staged PWM loop's grid, an injected
    row 1e-2 off from node 2 on; on the switched sampler the adjoint
    sideband rows (pnoise) missed the forward PAC by 1.1e-5 / 7.2e-5 at
    l = 0 / 1 (1e-15 opened at the manufacturing step).  Pinned: both to
    1e-10."""
    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric
    Tp = 1e-5
    cir = _pwm_loop(Tp)
    p = PSS(cir, method='trap', reltol=1e-10)
    with quiet(AccuracyWarning):
        ## (settled first since a landed ramp edge drops the order,
        ## 2026-09-28: from zeros trap's staged Newton stalls on this loop;
        ## `_staged_fallback` recovers it since, at ~9x the time)
        p.solve(period=Tp, timestep=Tp / 60, x0=np.zeros(cir.n - 1), maxiterations=100,
                tstab=20 * Tp)
    fp = p.factored_period()
    assert fp.is_plain and fp.open_at_x0
    m = cir.n - 1
    rng = np.random.default_rng(1)
    u, v = rng.standard_normal(m), rng.standard_normal(m)
    lhs = float(np.asarray(fp.matvec_transposed(u)) @ v)
    rhs = float(u @ np.asarray(fp.matvec(v)))
    assert abs(lhs - rhs) < 1e-10 * abs(rhs), (lhs, rhs)
    fclk = 100e3
    T = 1.0 / fclk
    fin = 7e3
    cir = SubCircuit()
    for nd in ('in', 'out', 'ck'):
        cir.add_node(nd)
    cir['Vin'] = VSin('in', gnd, vo=0.5, va=0.4, freq=fclk, phase=0.0, vac=1.0)
    cir['Vck'] = VSin('ck', gnd, vo=0.0, va=1.0, freq=fclk, phase=90.0)
    cir['S0'] = _SwitchHdl('in', 'out', 'ck', gnd, gon=1e-3, goff=1e-9, vth=0.0, vs=50e-3)
    cir['C0'] = C('out', gnd, c=100e-12)
    q = PSS(cir, method='trap', reltol=1e-10)
    with quiet(AccuracyWarning):
        q.solve(period=T, timestep=T / 200, x0=np.zeros(cir.n - 1), maxiterations=60,
                x0_unknown=True)
    assert q.factored_period().open_at_x0
    io_full = [str(n_) for n_ in cir.nodes].index('out')
    io = io_full if io_full < q.irefnode else io_full - 1
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), q.irefnode, circuit.numeric)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet():
        res = pac.solve(q, freqs=[fin])
        H = np.asarray(pac.adjoint_sideband_row(q, fin, io, sidebands=[0, 1])[0])
    fout = np.asarray(res.sweep_values, dtype=float)
    X = np.asarray(res.x)
    for li, l in enumerate((0, 1)):
        k = int(np.argmin(np.abs(fout - abs(fin + l * fclk))))
        x = complex(X[io_full, k])
        h = complex(H[li] @ u_ac)
        assert abs(x - h) < 1e-10 * abs(h), (l, x, h)


def test_a_twin_takes_the_analysis_reference_node():
    """The review of 2026-09-30 (F1): the monodromy twin of a trap/euler
    oscillator (and the staged fallback's one-stage solve) was built with
    `irefnode=None` -- ground -- while `solve` refuses a `refnode` other than
    the constructed one, so a PSS on another reference raised from `ppv()`.
    The diffusion constant is gauge-free: on the tank node as the reference
    it equals the ground-referenced one."""
    from pycircuit.circuit.tests.test_shooting_pss import _review_vdp
    cg, pg = _review_vdp()
    cv, pv = _review_vdp(ref='v')
    with quiet():
        c_v = PAC(cv, toolkit=circuit.numeric).diffusion_constant(pv)
        c_g = PAC(cg, toolkit=circuit.numeric).diffusion_constant(pg)
    assert pv.monodromy_twin().irefnode == pv.irefnode
    assert abs(c_v / c_g - 1.0) < 1e-9, (c_v, c_g)


def test_a_relative_sweep_offsets_from_the_carrier_of_the_map_it_solves_on():
    """The review of 2026-09-30 (F2): on a trap oscillator the small-signal
    surfaces solve on the monodromy twin's map, whose period differs from the
    run's by O(h^2) (3.4e-4 here, 5.4e-5 Hz at f0 = 0.159 Hz); a relative
    sweep read the RUN's carrier, so an offset of 1e-3 f0 sat 34 % off the
    twin's pole and `pnoise` read 2.3x off.  The offset is now from the
    carrier of the map the rows are solved on: `pnoise` at a relative
    offset IS `pnoise` at the twin's carrier plus that offset, and
    `PAC.solve`'s output frequencies are the twin's."""
    from pycircuit.circuit.tests.test_shooting_pss import _review_vdp
    cir, pss = _review_vdp(ac=True)
    tw = pss.monodromy_twin()
    f0r, f0t = 1.0 / float(pss.period), 1.0 / float(tw.period)
    df = 1e-3 * f0t
    ## (the poison is alive: the carriers are further apart than 10 % of df)
    assert abs(f0r - f0t) > 0.1 * df, (f0r, f0t, df)
    pac = PAC(cir, toolkit=circuit.numeric)
    with quiet(AccuracyWarning):
        S_rel, _u = pac.pnoise(pss, df, 0, maxsidebands=3, sweeptype='relative')
        S_abs, _u = pac.pnoise(pss, f0t + df, 0, maxsidebands=3,
                               sweeptype='absolute')
        res = pac.solve(pss, [df], sweeptype='relative')
    assert abs(S_rel / S_abs - 1.0) < 1e-12, (S_rel, S_abs)
    sv = np.asarray(res.sweep_values, dtype=float)
    assert float(np.min(np.abs(sv - (f0t + df)))) < 1e-9 * f0t, sv


def test_pac_solve_on_a_trap_staged_oscillator_reads_one_host():
    """The review of 2026-09-30 (F8): on a trap staged OSCILLATOR `PAC.solve`
    read the monodromy twin's map but the RUN's event columns, whose grids
    differ (the landed crossings) -- it crashed broadcasting 270 nodes
    against 236; `adjoint_sideband_row` read the run's columns beside the
    twin's map too.  Both take the twin whole now: the trap run's response
    IS its twin's, and meets an independent radau solve to O(h^2)."""
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    cir['ac'] = IS('c', gnd, i=0.0, iac=1e-6)
    seed, T = _relaxation_oscillator_seed(cir)
    runs = {}
    for method in ('trap', 'radau'):
        p = PSS(cir, method=method, reltol=1e-9)
        with quiet(AccuracyWarning):
            p.solve(period=T, timestep=T / 100, x0=seed, maxiterations=100)
        assert p.converged and p._event_columns is not None, method
        runs[method] = p
    tw = runs['trap'].monodromy_twin()
    assert tw is not runs['trap'] and tw._event_columns is not None
    io = [str(n) for n in cir.nodes].index('fb1')
    offs = np.array([1e-3, 0.1]) / float(tw.period)

    def response(p):
        pac = PAC(cir, toolkit=circuit.numeric)
        with quiet():
            res = pac.solve(p, offs, sweeptype='relative')
            row = pac.adjoint_sideband_row(p, 1.0 / float(tw.period) + offs[1],
                                           io, sidebands=[0, 1])[0]
        fo = np.asarray(res.sweep_values, float)
        X = np.asarray(res.x)
        f0p = 1.0 / float(p.monodromy_twin().period)
        return (np.array([X[io, int(np.argmin(np.abs(fo - (f0p + o))))] for o in offs]),
                np.asarray(row))

    r_run, h_run = response(runs['trap'])
    r_twin, h_twin = response(tw)
    r_rad, _h = response(runs['radau'])
    assert np.max(np.abs(r_run / r_twin - 1.0)) < 1e-12, r_run / r_twin - 1.0
    assert np.max(np.abs(h_run - h_twin)) <= 1e-12 * np.max(np.abs(h_twin))
    assert np.max(np.abs(r_run / r_rad - 1.0)) < 1e-4, r_run / r_rad - 1.0


def test_the_adjoint_transfer_row_tolerance_binds_on_the_iterative_path():
    """The review's X9 (2026-10-01): `adjoint_transfer_row(tol=)` (`recycle_tol`
    until batch 14) had
    no test -- and below `FLOQUET_DENSE_LIMIT` it CANNOT bind: the transposed
    solve is direct there (batch 4).  Forced onto the iterative path (the
    limit lowered on the PSS), it is the GMRES tolerance, monotone against
    the direct row on a 32-state ladder: 8.4e-3 / 3.4e-9 / 1.3e-15 at 1e-2 /
    1e-6 / 1e-12, in 5 / 9 / 11 mat-vecs (the default, KRYLOV_FACTOR *
    reltol: 1.4e-12)."""
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(30)
    pss = PSS(cir, method='gear', reltol=1e-9)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    f = 0.13e3
    direct, info = pac.adjoint_transfer_row(pss, f, 3)
    assert info['matvecs'] == 0
    pss.FLOQUET_DENSE_LIMIT = 0

    def err(tol):
        r, info = pac.adjoint_transfer_row(pss, f, 3, tol=tol)
        return (float(np.max(np.abs(r - direct)) / np.max(np.abs(direct))),
                info['matvecs'])
    loose, n_loose = err(1e-2)
    tight, n_tight = err(1e-12)
    assert loose > 1e-4 and tight < 1e-12 and n_loose < n_tight, \
        (loose, tight, n_loose, n_tight)


def test_pss_and_pac_run_under_the_jax_toolkit():
    """The review's X9 (2026-10-01): nothing ran the shooting stack under the
    JAX toolkit.  A driven RC at 10 points: the PSS waveform matches the
    numeric toolkit's to 1.2e-16 and the PAC sidebands exactly.  (Slow --
    per-call JAX dispatch in every element evaluation: ~12 s here against
    well under one numerically, and a diode mixer at 100 points did not
    finish in 15 minutes.)  `pnoise` too, plain and cyclostationary, equal
    to the numeric toolkit's: it failed until 2026-10-01 in `VS.CY`, which
    assigned into a JAX array in place (fails on the parent with "JAX arrays
    are immutable")."""
    pytest.importorskip('jax')
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit.toolkit import jaxtoolkit

    def run(tk):
        saved = _cm.default_toolkit
        _cm.default_toolkit = tk
        try:
            c = SubCircuit(toolkit=tk)
            c.add_node('1')
            c.add_node('2')
            c['vs'] = VSin('1', gnd, vac=1.0, va=1.0, freq=1e6,
                           noisePSD=1e-16, toolkit=tk)
            c['R'] = R('1', '2', r=1e3, noisy=True, toolkit=tk)
            c['C'] = C('2', gnd, c=1e-10, toolkit=tk)
            p = PSS(c, toolkit=tk, method='gear', reltol=1e-10)
            with quiet(AccuracyWarning):
                p.solve(period=1e-6, timestep=1e-6 / 10, maxiterations=20)
                pac = PAC(c, toolkit=tk)
                res = pac.solve(p, freqs=[0.13e6], sweeptype='absolute')
                S = [pac.pnoise(p, 0.13e6, '2', maxsidebands=3,
                                sweeptype='absolute', cyclostationary=cy)[0]
                     for cy in (False, True)]
            assert p.converged
            return (np.asarray(p.waveform[1], dtype=float), np.asarray(res.x),
                    np.asarray(S, dtype=float))
        finally:
            _cm.default_toolkit = saved
    Xn, Yn, Sn = run(circuit.numeric)
    Xj, Yj, Sj = run(jaxtoolkit)
    assert np.max(np.abs(Xj - Xn)) < 1e-12 * np.max(np.abs(Xn))
    assert np.max(np.abs(Yj - Yn)) < 1e-12 * np.max(np.abs(Yn))
    assert np.all(Sn > 0)
    assert np.max(np.abs(Sj - Sn)) < 1e-12 * np.max(Sn)
