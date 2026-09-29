"""Shooting tests: shooting events.  Split out of test_analysis_shooting.py on
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
from pycircuit.circuit.tests._shooting_elements import _SwitchHdl
from pycircuit.circuit.tests._shooting_fixtures import (_comparator_relaxation_oscillator,
    _e5_circuit,
    _pulse_clocked_sampler,
    _pwm_loop,
    _relaxation_oscillator_seed)


def test_event_grid_lands_the_period_on_its_event_times():
    """A6/B7b: `PSS.event_grid` puts the circuit's event times ON grid points.

    `Transient` breaks its steps at `cir.next_event`; the PSS traversal does not,
    so a pulse edge inside a step is integrated straight through. Landing the
    edges costs 3-4 points out of 40 and buys up to 27x accuracy.

    ⚠ THE SNAP IS WHY THERE ARE NO SLIVERS. An event near an existing point MOVES
    that point onto it rather than inserting a second one beside it; inserting
    unconditionally is how a merge acquires arbitrarily small steps, which is the
    same lesson B7c's separation rule encoded (`refine_grid`, deleted
    2026-09-21 -- `lte_grid` from a converged state is the repair path).
    Check 3 pins it.

    ⚠ AND THIS IS NOT SALTATION, which was measured and falsified twice for this
    codebase -- a switched conductance, a discontinuous injection and an
    `Idtmod` wrap all give a monodromy-vs-FD gap falling at 2.00x per doubling,
    i.e. O(h). Each step already uses its own `Jf`/`C`, describing whichever side
    of the switch it is on. The defect was only ever that the grid could not
    BREAK at the event.

    ⚠ TIME-DRIVEN EVENTS ONLY. A state-dependent reset cannot be walked out in
    advance (`Idtmod.next_event` is a linear prediction from the last accepted
    point and is `inf` before a traversal starts), and its wrap time MOVES as the
    Newton iterates. That case needs the event time to become a Newton unknown
    and is deliberately out of scope here.
    """
    import warnings
    from pycircuit.circuit.elements import VPulse
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    N = 40

    def cir(td):
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=td, tr=T / 500,
                         tf=T / 500, pw=T / 2, per=T)
        c['R'] = R(1, 2, r=1e3)
        c['C'] = C(2, gnd, c=1e-9)
        return c

    def solve(fr, cr):
        p = PSS(cr, method='gear', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ## the grid under test is the one passed in -- since 2026-09-20
            ## `break_events` defaults ON for gear too, and would land the
            ## events on the "uniform" baseline as well (measured: identical)
            p.solve(period=T, grid=list(fr), maxiterations=40,
                    break_events=False)
        return np.asarray(p._period_state[1], float).ravel()

    td = 0.5 / N * T                      # edges deliberately between points
    uni = list(np.diff(np.linspace(0.0, 1.0, N + 1)))
    ref = solve(list(np.diff(np.linspace(0.0, 1.0, 4001))), cir(td))
    den = max(np.max(np.abs(ref)), 1e-30)
    err_uni = np.max(np.abs(solve(uni, cir(td)) - ref)) / den

    p = PSS(cir(td), method='gear', reltol=1e-10)
    eg = p.event_grid(T, npts=N)
    err_ev = np.max(np.abs(solve(eg, cir(td)) - ref)) / den

    ## 1. the events are actually ON the grid
    pts = np.concatenate(([0.0], np.cumsum(np.asarray(eg))))
    assert p.event_times, 'no events were found on a VPulse-driven circuit'
    for f in p.event_times:
        assert np.min(np.abs(pts - f)) < 1e-12, \
            'event at %.9f of the period is not on a grid point' % f

    ## 2. and landing them is worth doing
    assert err_ev < err_uni / 3.0, \
        'landing the edges gave %.3e against the uniform grid\'s %.3e, which ' \
        'is less than the 3x this costs points for' % (err_ev, err_uni)

    ## 3. NO SLIVERS -- the snap must keep every step a real fraction of one
    h_uni = 1.0 / N
    assert min(eg) > 0.02 * h_uni, \
        'the smallest step is %.3e, i.e. %.1f%% of a uniform step -- the snap ' \
        'is not preventing slivers' % (min(eg), 100.0 * min(eg) / h_uni)

    ## 4. a circuit with no events is left exactly alone
    c2 = SubCircuit()
    c2['vs'] = VSin(1, gnd, va=1.0, freq=1.0 / T)
    c2['R'] = R(1, 2, r=1e3)
    c2['C'] = C(2, gnd, c=1e-9)
    p2 = PSS(c2, method='gear', reltol=1e-10)
    base = list(np.diff(np.linspace(0.0, 1.0, 11)))
    got = p2.event_grid(T, grid=base)
    if not p2.event_times:
        assert np.allclose(got, base, rtol=0, atol=1e-15), \
            'a circuit with no events had its grid altered'


def test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_not_on_the_grid(monkeypatch):
    """Where an `Idtmod` wrap falls RELATIVE TO THE GRID does not matter; where
    the FINAL state falls relative to the fold is the whole story.

    ⚠ THIS TEST EXISTS TO RECORD A FALSIFICATION.  The roadmap (A6) recorded
    that the period map is "genuinely discontinuous" when the reset lands ON a
    grid point -- `|dphi|` a constant 8.2e-3 independent of eps at npts=1200
    against a smooth 1.732 at npts=600 -- and reasoned that "an infinitesimal
    change flips which step the reset lands in and a whole modulus propagates".
    On that reading the remedy was to make the event time an unknown the Newton
    solves for, inside the transient engine, and to split each monodromy step
    into two composed sub-step maps.  That is expensive surgery, so it was
    priced with a falsifier first: does landing the wrap EXACTLY remove the
    jump?

    It does not, and the premise does not survive either:

    * The jump is **grid-independent**.  Measured at npts 250/500/600/1000/
      1200/2000/2400 with the wrap at node 86.25 (off), 345.00 (exactly on),
      and five others: `||dphi||/eps` is 1.414214e+09 at EVERY ONE, to all
      printed digits.  A quantity that does not move when the grid moves is not
      a grid artefact.
    * Adding the exact wrap times to the grid leaves it at 1.414214e+09.
    * A dense scan of 60 base offsets across (0.001, 0.499), at npts 600 and
      1200 -- every one of which places the wrap on an integer node -- found
      **zero** jumps and **zero** disagreements between the two grids.
    * The one base that jumps is the one where the final state sits on the
      fold: `phi_idt(T)` steps from -0.000000000000 to -0.999999999999 as the
      offset crosses zero.  That is the modulus, applied at t = T.

    So the discontinuity set of the period map is `{x0 : phi(x0) lands on a
    fold boundary}` -- a measure-zero set fixed by the OUTPUT map, which no
    refinement of the time grid can move.  Event localisation is the wrong
    instrument for it; the right one is to stop differencing across the fold
    (carry the unfolded phase, or take the shooting residual modulo the
    modulus).  A6's time-driven half stands as built -- `event_grid` is
    measured to help real sources -- but its state-dependent half is NOT the
    item this record claimed, and the traversal surgery is not justified by it.

    ⚠ Caveat kept deliberately: the recorded fixture was not reproduced
    exactly.  Its smooth column reads 1.732 = sqrt(3) where this one reads
    sqrt(2), so the recorded circuit had a third responding coordinate that
    this one does not.  What is asserted here is what THIS fixture measures.
    """
    from copy import copy
    circuit.default_toolkit = circuit.numeric

    T, IC, RATE = 1e-3, 0.31, 2.0

    def build():
        c = SubCircuit()
        c['vin'] = VS('in', gnd, v=RATE / T)
        c['X'] = Idtmod('in', gnd, 'o', gnd, modulus=1.0, ic=IC)
        c['Ro'] = R('o', gnd, r=1e6)
        return c

    n = build().n
    idt = 1                       # reduced coord of 'X.idt_node'

    def phi(x0r, times):
        p = PSS(build(), method='gear', reltol=1e-11)
        p._tran = p._new_transient(p._integrator_for('gear'))
        p._want_dfdh = False
        p._want_lte = False
        p._begin_period(np.asarray(x0r, dtype=float))
        x = copy(np.asarray(x0r, dtype=float))
        hs = np.diff(times)
        for j, t in enumerate(times[1:]):
            x = copy(p.solve_timestep(x, t, hs[j]))
        return np.asarray(x, dtype=float).ravel()

    ## ⚠ THE STAGE PREDICTOR IS PINNED OFF, and the reason is this fixture's
    ## own subject.  The jump measured below IS a fold -- the map is
    ## DISCONTINUOUS there -- so which side of it a given grid's endpoint lands
    ## on is decided at the last ulp of the trajectory.  A predictor changes the
    ## trajectory by about that much (measured, 3e-14 on a converged orbit) and
    ## one grid of the four then lands the other side: the jump reads 2.000e9
    ## instead of 2.236e9, both of them jumps, neither of them a grid artefact.
    ## The QUALITATIVE claim -- every grid jumps -- is re-checked with the
    ## predictor on at the end, which is the part that is not knife-edge.  ⚠ The
    ## history of this fixture is that the FIRST one sat exactly ON the fold and
    ## would have confirmed the opposite story; it is fragile by construction.
    from pycircuit.circuit.transient import Transient
    monkeypatch.setattr(Transient, 'stage_predictor', 'off')

    d = np.zeros(n - 1)
    d[idt] = 1.0
    eps = 1e-9

    ## The perturbation direction must actually drive the map, or every number
    ## below is a zero-vs-zero pass.
    g = np.linspace(0, T, 1201)
    base = phi(0.20 * d, g)
    assert np.linalg.norm(phi(0.20 * d + eps * d, g) - base) / eps > 1.0, \
        'the idt state does not propagate -- the fixture proves nothing'

    ## (1) ON a grid node is not special.  Every offset here puts the first
    ## wrap at t/T = (1 - IC - b)/RATE, i.e. exactly on a node of a 1200-point
    ## grid, and every one is smooth.
    for b in (0.05, 0.10, 0.20, 0.35, 0.45):
        r = np.linalg.norm(phi(b * d + eps * d, g) - phi(b * d, g)) / eps
        assert abs((1 - IC - b) / RATE * 1200 - round((1 - IC - b) / RATE * 1200)) < 1e-9, \
            'b=%g does not place the wrap on a node -- the premise is untested' % b
        assert r < 10.0, \
            'wrap on node %d: expected a smooth map, got ||dphi||/eps = %.4e' \
            % (round((1 - IC - b) / RATE * 1200), r)

    ## (2) The jump is at the ENDPOINT fold, and it does not care about the grid.
    ratios = []
    for npts in (250, 600, 1000, 1200):
        gg = np.linspace(0, T, npts + 1)
        ratios.append(np.linalg.norm(phi(eps * d, gg) - phi(np.zeros(n - 1), gg)) / eps)
    assert min(ratios) > 1e8, \
        'the endpoint fold should jump on every grid, got %s' % ratios
    assert max(ratios) - min(ratios) < 1e-3 * max(ratios), \
        'the jump moved with the grid (%s) -- it WOULD then be a grid artefact ' \
        'and A6\'s recorded reading would stand' % ratios


    ## and with the predictor on, the qualitative claim: still a jump on every
    ## grid, still ~1e9 (the modulus over eps), just not comparable digit for
    ## digit across grids at a discontinuity
    monkeypatch.setattr(Transient, 'stage_predictor', 'on')
    on_ratios = []
    for npts in (250, 600, 1000, 1200):
        gg = np.linspace(0, T, npts + 1)
        on_ratios.append(np.linalg.norm(phi(eps * d, gg)
                                        - phi(np.zeros(n - 1), gg)) / eps)
    assert min(on_ratios) > 1e8, on_ratios
    assert max(on_ratios) < 1e10, on_ratios
    ## (3) Landing the wrap exactly -- the falsifier -- does not help.
    wraps = [(k - IC) / RATE * T for k in (1, 2) if 0.0 < (k - IC) / RATE < 1.0]
    ge = np.unique(np.r_[np.linspace(0, T, 1201), wraps])
    r_exact = np.linalg.norm(phi(eps * d, ge) - phi(np.zeros(n - 1), ge)) / eps
    assert r_exact > 1e8, \
        'exact event placement removed the jump (%.4e) -- if this ever fires, ' \
        'the traversal surgery IS justified and this record must be reopened' % r_exact

    ## (4) The mechanism, asserted rather than described: a whole modulus at t=T.
    lo = phi(-1e-12 * d, g)[idt]
    hi = phi(+1e-12 * d, g)[idt]
    assert abs(abs(hi - lo) - 1.0) < 1e-6, \
        'expected exactly one modulus across the fold, got %.6e' % abs(hi - lo)


def test_the_shooting_residual_folds_a_periodic_state_and_leaves_everything_else_alone():
    """`x_0 - phi(x_0) == 0` is the WRONG condition on a state that is only
    defined up to `n*modulus`, and the circuit already says which states those
    are.

    `Idtmod` declares its integral row through `periodic_states()` -- the same
    declaration `Transient` uses for its gauge shift -- but `shooting.py` did
    not consume it, exactly as it did not consume `next_event`.  So shooting
    asked a folding row for the SAME REPRESENTATIVE rather than the same state:
    an orbit that closes after advancing one modulus had no root at all, and
    near the fold the raw difference jumped a whole modulus while the state
    moved infinitesimally.

    This is the residual-side repair, and it is the RIGHT instrument because
    the defect was measured to be in the output map, not the time grid --
    see `test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_not_on_the_grid`,
    where the jump is identical across seven grids.

    What is asserted, and what is deliberately NOT:

    * the gauge is collected and mapped to REDUCED rows, with the declared
      window's offset (the state row's fold needs only the modulus -- a
      DIFFERENCE is defined up to n*modulus whatever window each state sits
      in -- but `_wrap_jump` reads the window's edges, where the output wraps);
    * a circuit with no folding state gets an empty gauge and a bit-identical
      solution, so this cannot perturb the rest of the suite;
    * on a folding orbit the fold FIRES with a full-modulus correction and
      turns a LinAlgError into a converged, genuine orbit;
    * ⚠ it does NOT fix a free-running phase's singular `I - M`.  A DC-driven
      integrator has monodromy eigenvalue exactly 1 on its phase row, so at a
      rate of a whole number of moduli EVERY x_0 is a solution and the
      fixed-period Jacobian is singular *correctly* -- the problem is
      underdetermined, and locking it needs feedback (a PLL), not a residual
      change.  Asserted here so the fold is not credited with more than it does.
    """
    circuit.default_toolkit = circuit.numeric
    T = 1e-3

    def folding(rate):
        c = SubCircuit()
        c['vin'] = VS('in', gnd, v=rate / T)
        c['X'] = Idtmod('in', gnd, 'o', gnd, modulus=1.0, ic=0.31)
        c['Ro'] = R('o', gnd, r=1e6)
        return c

    ## (1) The gauge: global row -> reduced row, with the declared window
    ## (`Idtmod` declares ``-(offset + modulus)``, here -1).
    c = folding(0.5)
    declared = c.periodic_states()
    assert len(declared) == 1 and declared[0][1] == 1.0, \
        'the fixture stopped declaring a periodic state: %r' % (declared,)
    grow = int(declared[0][0])
    p = PSS(folding(0.5))
    gauge = p._collect_periodic_fold()
    iref = p.irefnode
    assert gauge == [(grow if grow < iref else grow - 1, 1.0, -1.0)], \
        'gauge %r does not map global row %d through irefnode %d' \
        % (gauge, grow, iref)

    ## (2) The fold itself: into [-m/2, m/2), identity well inside it.
    p._periodic_fold = gauge
    r = gauge[0][0]
    probe = np.zeros(p.cir.n - 1)
    probe[r] = 1.0
    assert abs(p._fold_periodic(probe)[r]) < 1e-12, \
        'a full modulus must fold to zero -- that is the whole point'
    probe[r] = 0.25
    assert abs(p._fold_periodic(probe)[r] - 0.25) < 1e-15, \
        'the fold must be the identity away from the boundary'
    probe[r] = -1.75
    assert abs(p._fold_periodic(probe)[r] - 0.25) < 1e-12, \
        'the fold must reach the nearest representative, not just one modulus'

    ## (3) A circuit with no folding state is untouched, bit for bit.
    def plain():
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, va=2.0, freq=1e6, phase=20)
        c['R'] = R(1, 2, r=1e4)
        c['D'] = Diode(2, gnd)
        c['C'] = C(2, gnd, c=1e-12)
        return c
    assert plain().periodic_states() == []
    assert PSS(plain())._collect_periodic_fold() == []
    got, conv = [], []
    for fold in (False, True):
        q = PSS(plain(), method='gear')
        if not fold:
            q._collect_periodic_fold = lambda: []
        res = q.solve(period=1e-6, timestep=1e-6 / 400)
        got.append(np.asarray(q._period_state[1], dtype=float).ravel())
        conv.append(bool(q.converged))
    ## ⚠ Comparing two solves that BOTH failed would make this vacuous -- a
    ## pair of identical non-answers is still identical.
    assert all(conv), 'the invariance check needs converged solves, got %r' % conv
    assert np.array_equal(got[0], got[1]), \
        'the fold perturbed a circuit that declares no periodic state ' \
        '(||d|| = %.3e)' % np.linalg.norm(got[1] - got[0])

    ## (4) On a folding orbit it fires, and what it converges to is real.
    ## rate = 0.5 moduli per seed period: the orbit closes only after TWO,
    ## which is precisely the closure the unfolded residual cannot express.
    import warnings as _w
    from numpy.linalg import LinAlgError

    def attempt(fold):
        q = PSS(folding(0.5), method='trap')
        if not fold:
            q._collect_periodic_fold = lambda: []
        hits = {'n': 0, 'max': 0.0}
        orig = q._fold_periodic

        def spy(F):
            G = orig(F)
            d = np.max(np.abs(np.asarray(F, dtype=float).ravel()
                              - np.asarray(G, dtype=float).ravel()))
            if d > 1e-9:
                hits['n'] += 1
                hits['max'] = max(hits['max'], d)
            return G
        q._fold_periodic = spy
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=T, timestep=T / 300)
        return q, hits

    with pytest.raises(LinAlgError):
        attempt(False)

    q, hits = attempt(True)
    assert hits['n'] >= 1 and abs(hits['max'] - 1.0) < 1e-6, \
        'the fold did not fire a full-modulus correction (%r) -- then it is ' \
        'not what rescued this solve' % (hits,)

    ## ⚠ Converged is not solved.  Re-traverse from the solution and require
    ## the OPENED state (not the raw unknown, whose algebraic entries are
    ## free) to return to itself on EVERY row.
    solved, x_in, xm1, times, hs, Tp, x0u = q._period_state
    wk = q._walk('plain', np.asarray(x_in, dtype=float).ravel(), times, hs,
                 T=Tp, open_at_x0=x0u)
    x0, x_end = wk.x0, wk.x_end
    gap = np.max(np.abs(np.asarray(x0, dtype=float).ravel()
                        - np.asarray(x_end, dtype=float).ravel()))
    assert gap < 1e-10, \
        'the folded residual converged to a point that is NOT an orbit ' \
        '(max row gap %.3e) -- a fold that admits spurious roots is worse ' \
        'than no fold' % gap
    assert abs(Tp / T - 2.0) < 1e-6, \
        'expected the orbit to close after two seed periods, got T/T0 = %.6f' \
        % (Tp / T)

    ## (5) The limit of the claim, asserted rather than asserted-about: the
    ## phase row is a MARGINAL mode.  Perturbing it moves the endpoint by
    ## exactly the same amount, so the period map's derivative along it is 1
    ## and `I - M` is singular there by construction.  That is why a driven
    ## fixed-period solve on a free-running integrator is underdetermined, and
    ## no residual fold can or should change it.
    from copy import copy as _copy
    probe = PSS(folding(0.5), method='gear')
    probe._tran = probe._new_transient(probe._integrator_for('gear'))
    probe._want_dfdh = False
    probe._want_lte = False
    tt = np.linspace(0, T, 401)

    def endpoint(v):
        probe._begin_period(np.asarray(v, dtype=float))
        x = _copy(np.asarray(v, dtype=float))
        steps = np.diff(tt)
        for j, t in enumerate(tt[1:]):
            x = _copy(probe.solve_timestep(x, t, steps[j]))
        return np.asarray(x, dtype=float).ravel()

    base = np.zeros(probe.cir.n - 1)
    base[r] = 0.2                     # off the fold boundary
    d = np.zeros(probe.cir.n - 1)
    d[r] = 1.0
    eps = 1e-7
    slope = (endpoint(base + eps * d) - endpoint(base))[r] / eps
    assert abs(slope - 1.0) < 1e-5, \
        'the phase row should be marginal (dx_end/dx_0 == 1 on it), got ' \
        '%.9f -- if this ever stops being 1 the underdetermination argument ' \
        'above needs rewriting' % slope


def test_an_orbit_that_starts_on_an_output_wrap_closes_across_the_jump():
    """The state fold leaves an idtmod's OUTPUTS discontinuous: when the
    orbit wraps exactly at ``t = 0``, `z_0` sits just after the wrap and
    `z_end` just before it, and every algebraic quantity the wrapped output
    feeds differs by the whole jump.  A free-running `VcoHdl` pinned at a zero
    crossing of `sin(2 pi phase)` starts exactly there, and radau failed on
    it at every load (2026-09-23): the Newton iterates straddled the wrap at
    rounding level and the jump flipped in and out of the residual.
    `_wrap_jump` subtracts the orbit's own jump, decided by the STATES.

    The loads matter: a fold of the phase node by its modulus closes the
    unloaded node and nothing else -- a load resistor's current jumps by
    ``modulus/R`` and a gain-2 detector's node by ``2*modulus`` (measured:
    both still failed with the node folded).  And the fold must be decided
    by the state, not the residual: rounding the residual by the modulus
    accepted a seed one modulus off its constraint (`ph(0) = 1.1` against a
    phase of 0.1) as a converged orbit under radau, gear and trap.
    """
    circuit.default_toolkit = circuit.numeric
    from pycircuit.circuit.elements_hdl import VcoHdl

    def vco(load):
        c = SubCircuit()
        for nd in ('vco', 'ph', 'ctl') + (('pd',) if load == 'E' else ()):
            c.add_node(nd)
        c['X1'] = VcoHdl('ctl', gnd, 'vco', gnd, 'ph', f0=1e6, kvco=0.0,
                         va=1.0, modulus=1.0)
        c['Rc'] = R('ctl', gnd, r=1e3)
        c['Rl'] = R('vco', gnd, r=1e3)
        if load == 'R':
            c['Rph'] = R('ph', gnd, r=1e3)
        else:
            c['E1'] = VCVS('ph', gnd, 'pd', gnd, g=2.0)
            c['Rpd'] = R('pd', gnd, r=1e3)
        return c

    def solve(load, x0):
        c = vco(load)
        names = [str(nd) for nd in c.nodes if str(nd) != 'gnd!']
        p = PSS(c, method='radau', reltol=1e-9)
        p.solve(period=1.03e-6, timestep=1.03e-6 / 100, x0=x0(c, names),
                maxiterations=30)
        return p, names, np.asarray(p.waveform[1], dtype=float)

    ## (1) Seeded consistently ON the wrap (every entry zero): the orbit
    ## starts there, and closes.  ⚠ Whether the iterates straddle is
    ## rounding luck: at 40 points the code before the fix converged here, at
    ## 100 (this grid) it failed on both loads.
    for load in ('R', 'E'):
        p, names, xs = solve(load, lambda c, names: np.zeros(c.n - 1))
        ip, ist = names.index('ph'), names.index('X1._state0')
        assert p.converged, '%s: the orbit starting on the wrap did not close' % load
        assert abs(p.period - 1e-6) < 1e-15, (load, p.period)
        ## ⚠ Not vacuous: the solved orbit really does start on the wrap.
        assert abs(xs[ip, 0]) < 1e-12, \
            '%s: ph(0) = %.3e -- the orbit no longer starts on the wrap, so ' \
            'this no longer tests the jump' % (load, xs[ip, 0])
        d = (xs[ip] - np.mod(xs[ist], 1.0) + 0.5) % 1.0 - 0.5
        assert np.max(np.abs(d)) < 1e-12, \
            '%s: ph is not the wrap of its state (%.3e)' % (load, np.max(np.abs(d)))

    ## (2) The jump itself, read at a hand-built straddle: exactly the output
    ## map's jump on every row it feeds, nothing on the state's own row, and
    ## nothing at all when the two states sit on the same branch.
    ## (the waveform carries the reference row; the residual does not -- the
    ## nodes named here all precede it)
    k = {nm: j for j, nm in enumerate(names)}
    z0 = np.delete(xs[:, 0], p.irefnode)
    ze = z0.copy()
    z0[k['X1._state0']], ze[k['X1._state0']] = 1e-9, 1.0 - 1e-9
    jump = p._wrap_jump(z0, ze, (0.0, p.period))
    assert abs(jump[k['ph']] + 1.0) < 1e-12, jump
    assert abs(jump[k['pd']] + 2.0) < 1e-12, jump
    ## (`vco = sin(2 pi ph)` is continuous: what is left is the 2x64-ulp
    ## straddle the jump is read across)
    assert abs(jump[k['X1._state0']]) < 1e-15 and abs(jump[k['vco']]) < 1e-12, jump
    z0[k['X1._state0']], ze[k['X1._state0']] = 0.3, 1.3
    assert not np.any(p._wrap_jump(z0, ze, (0.0, p.period)))

    ## (3) No false root: a seed one modulus off its own constraint is
    ## repaired, not accepted.
    def off_by_one(c, names):
        x0 = np.zeros(c.n - 1)
        x0[names.index('X1._state0')] = 0.1
        x0[names.index('vco')] = np.sin(0.2 * np.pi)
        x0[names.index('ph')] = 1.1
        return x0
    p, names, xs = solve('R', off_by_one)
    assert p.converged
    assert abs(xs[names.index('ph'), 0] - 0.1) < 1e-9, xs[names.index('ph'), 0]


def test_event_breaking_defaults_on_for_every_method_and_a_jump_keeps_both_ramp_ends():
    """The default is ON for every method, and `event_grid` never snaps an
    event onto an event.  Both are measured on a ladder (2026-09-20).

    ⚠ THIS USED TO ASSERT GEAR DEFAULTS OFF, "measured to lose", on ONE step
    count (gear uniform 8.23e-3, + events 1.29e-2, with a jittered control
    showing non-uniformity itself costs gear).  The ladder overturned the
    conclusion, not the control: non-uniformity does cost a multistep method,
    and the alignment buys more.  Pulsed RC, analytic reference, max error
    over the ACTUAL grid nodes, edges ramped over 0.02 T::

        N       gear uniform  gear + events   trap uniform  trap + events
        50      1.39e-02      3.16e-03        7.78e-03      1.68e-04
        400     3.72e-04      2.44e-04        1.30e-04      5.81e-06
        800     9.52e-05      6.31e-05        3.25e-05      1.45e-06
        1600    2.41e-05      1.61e-05        8.13e-06      3.63e-07

    gear + events wins at 5 of 6 N and is second order (3.88x, 3.92x per
    doubling).  Its gain (1.5x) is small next to trap's (22x) because BDF2
    takes an O(h^2 [x'']) hit at the ONE step after each corner, where its
    history straddles the jump in x'' (5.5e-7 -> 1.7e-4 across that step,
    30x trap's) -- a constant, not an order.

    ⚠ AT A TRUE JUMP (tr = 0, clamped to a 1e-18 ramp) EVERY METHOD WAS
    FIRST ORDER WITH THE EDGES LANDED, and the cause was the grid, not a
    method: the snap landed the ramp's first event and then OVERWROTE that
    node with its second, so the edge sat on the post-jump side and the
    step arriving there integrated its whole length with the post-jump
    source.  The error was injected AT the edge node and decayed over tau,
    each method's endpoint weight times h dU/tau (radau 0.111 x 8.3e-3 =
    9.3e-4 predicted, 8.5e-4 measured).  Keeping BOTH ramp ends as nodes --
    a 1e-18 step -- fixes all three, gear included, with no restart::

        N      gear collapsed  gear both ends   trap both ends  radau both ends
        200    5.73e-03        6.43e-05         6.27e-06        5.0e-08
        800    2.55e-03        1.01e-05         6.34e-07        2.0e-07
        1600   1.30e-03        3.42e-06         2.30e-06        5.3e-07

    Variable-step BDF2 with h_n/h_{n-1} -> inf over a consistent tiny step
    degenerates to the trapezoidal rule (the factor w/2 multiplies a
    difference that is O(1/w)), which is also why the zero-stability
    warning must not fire on an ISOLATED up-step.  The ~1e-6 floors at
    800-1600 are the ramp node's time rounding (1e-24 in t is 1e-6 in u on
    a 1e-18 ramp), not the solver.

    ⚠ The reference is the RC's closed form per linear segment written with
    `expm1`: the naive `a - b tau` with b = 1/tr = 1e18 cancelled to a flat
    9.5e-5 on every method, and that flat floor was the tell.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    RR, CC = 1e3, 3e-10
    TAU = RR * CC
    TD, PW = 0.0125 * T, 0.4 * T

    def pulsed(tr):
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=TD, tr=tr, tf=tr, pw=PW, per=T)
        c['R'] = R(1, 2, r=RR)
        c['C'] = C(2, gnd, c=CC)
        return c

    def exact(tr, ts):
        tr = max(tr, 1e-18)
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

    def err(tr, method, N, be):
        c = pulsed(tr)
        pss = PSS(c, method=method, reltol=1e-10)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            pss.solve(period=T, timestep=T / N, break_events=be, maxiterations=40)
        assert pss.converged
        ts = np.asarray(pss.waveform[0], float).ravel()
        xs = np.asarray(pss.waveform[1], float)
        row = [str(n) for n in c.nodes].index('2')
        ratio_warn = [w for w in rec if 'steps up by' in str(w.message)]
        return float(np.abs(xs[row] - exact(tr, ts)).max()), pss, ratio_warn

    ## (1) ON for every method; an explicit value honoured both ways.
    for meth in ('euler', 'trap', 'gear', 'radau', 'trbdf2'):
        assert PSS(pulsed(0.02 * T), method=meth)._resolve_break_events(None) is True, meth
    assert PSS(pulsed(0.02 * T), method='gear')._resolve_break_events(False) is False
    assert PSS(pulsed(0.02 * T), method='trap')._resolve_break_events(True) is True

    ## (2) Ramped edges: gear + events second order.
    ## ⚠ NO LONGER AHEAD OF UNIFORM ON THIS FIXTURE (2026-09-28): a landed
    ## edge now drops the order (the stepping loop's rule, `solve_timestep`),
    ## which on a SMOOTH state costs more than the landing buys -- gear with
    ## events 6.60e-4 against uniform 3.72e-4 at 400 (1.69e-4 / 9.52e-5 at
    ## 800), trap 4.61e-4 against 1.30e-4 -- while on a STIFF one it is the
    ## only grid that does not ring (uniform trap 0.83 of the current there);
    ## `test_the_landed_edges_drop_the_order_so_a_stiff_state_does_not_ring`.
    e_e800, _, _ = err(0.02 * T, 'gear', 800, True)
    e_e1600, g1600, wr = err(0.02 * T, 'gear', 1600, True)
    assert g1600.break_events is True and len(g1600.event_times) == 4
    assert e_e1600 < 5e-5, e_e1600
    assert 3.3 < e_e800 / e_e1600 < 4.5, (e_e800, e_e1600)
    assert not wr, 'an isolated event insertion must not raise the ratio warning'

    ## (3) A true jump: both ramp ends are nodes, and every method is back
    ## to its order.  The collapsed grid gave >= 1.3e-3 for all three.
    fr = PSS(pulsed(0.0)).event_grid(T, npts=200)
    assert min(fr) < 1e-9, 'the 1e-18 ramp step must survive: min step %.3e' % min(fr)
    assert len(fr) >= 204, len(fr)
    e_g200, _, _ = err(0.0, 'gear', 200, True)
    e_g800, g800, wr = err(0.0, 'gear', 800, True)
    assert e_g800 < 3e-5 and e_g200 / e_g800 > 3.0, (e_g200, e_g800)
    assert not wr, 'a 1e-18 event ramp is an isolated up-step: no ratio warning'
    e_t800, _, _ = err(0.0, 'trap', 800, True)
    ## (7.5e-6 since trap drops the order at both ramp ends, 2026-09-28:
    ## 6.4e-7 without; gear's growth guard already dropped there)
    assert e_t800 < 1e-5, e_t800
    e_r400, _, _ = err(0.0, 'radau', 400, True)
    assert e_r400 < 1e-7, e_r400

    ## (4) The ratio warning still fires where up-steps REPEAT.
    N = 200
    alt = np.tile([3.0, 1.0], N // 2)
    c = pulsed(0.02 * T)
    pss = PSS(c, method='gear', reltol=1e-8)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        pss.solve(period=T, timestep=T / N, break_events=False,
                  grid=alt / alt.sum(), maxiterations=40)
    assert any('steps up by' in str(w.message) and 'REPEAT' in str(w.message)
               for w in rec), 'an alternating 3:1 grid must still warn'

    ## (5) ⚠ A circuit whose sources declare no discontinuity must be
    ## BIT-IDENTICAL either way, or this default silently moves every solve in
    ## the suite.  `event_grid` rebuilds a uniform grid from `linspace` even
    ## when it finds nothing, and that differs from `_period_grid`'s in the
    ## last bit -- which is why `solve` only swaps the grid when an event
    ## actually exists.
    def smooth():
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, va=2.0, freq=1e6, phase=20)
        c['R'] = R(1, 2, r=1e4)
        c['D'] = Diode(2, gnd)
        c['C'] = C(2, gnd, c=1e-12)
        return c
    assert PSS(smooth()).event_grid(T, npts=40) is not None
    same = []
    for be in (False, True):
        q = PSS(smooth(), method='trap')
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=T, timestep=T / 200, break_events=be)
        same.append(np.asarray(q._period_state[1], dtype=float).ravel())
        assert q.event_times == [], \
            'the smooth fixture grew an event: %r' % (q.event_times,)
    assert np.array_equal(same[0], same[1]), \
        'break_events perturbed an event-free circuit (||d|| = %.3e)' \
        % np.linalg.norm(same[1] - same[0])


def test_algebraic_conditioning_resolves_a_crossing_rank_C_cannot_see():
    """The case a RANK test cannot answer, and the width of the window where
    this one cannot either.

    A binary "is `N^T G N` singular" test is blind here twice over: `rank C`
    never moves, and with a ONE-dimensional null space the usual
    `s.min()/s.max()` guard is identically 1, so it can only fire when the
    block is EXACTLY zero.  `algebraic_conditioning` returns the number
    instead, and tracks `G_22` down four decades.

    ⚠⚠ IT DOES HAVE A FALSE WINDOW -- every numerical index test does -- AND
    THE POINT IS WHERE IT IS.  An absolute rank test on these blocks has a
    window of width `~ tol * ||G|| / ||C||`, which WIDENS AS `1/||C||`: at
    picofarads it swallows perfectly healthy operating points, silently.
    This one's window is set by `||G||` ALONE, at about `1e-4 * ||G||`, and is
    invariant to the capacitance unit across twelve decades.  That invariance
    is the property under test; the width itself is just a number to know.
    """
    ## tracks the real value down four decades
    for g22 in (2e-3, 1e-3, 2e-4, 2e-5, 2e-6):
        sigma, info = algebraic_conditioning(_e5_circuit(g22))
        assert info['verdict'] == 'well-conditioned', (g22, info['verdict'])
        assert abs(sigma - g22) <= (info['spread'] - 1.0) * g22, (sigma, g22)
        ## and rank C never moved -- this is not a topology change
        assert topological_index(_e5_circuit(g22))[0] == 1

    ## exactly at the crossing it is genuinely singular
    sigma, info = algebraic_conditioning(_e5_circuit(0.0))
    assert info['verdict'] == 'singular', info['verdict']
    assert sigma == 0.0

    ## ⚠ THE INVARIANCE IS THE GATE.  Find the last resolved `G_22` at each
    ## capacitance scale; if the window were `1/||C||` this walks by twelve
    ## decades.  It does not move at all.
    limits = []
    for cscale in (1e6, 1e3, 1.0, 1e-3, 1e-6):
        last = None
        for e in range(2, 13):
            g22 = 2.0 * 10.0 ** (-e)
            _s, info = algebraic_conditioning(_e5_circuit(g22, cscale))
            if info['verdict'] != 'well-conditioned':
                break
            last = g22
        limits.append(last)
    assert len(set(limits)) == 1, limits

    ## ⚠⚠ AND THE WIDTH ITSELF IS A PARAMETER ARTEFACT, SO PIN THE CAUSE AND
    ## NOT THE NUMBER.  A first version of this gate asserted `== 2e-7` flat,
    ## which reads as a property of the problem and is not one: it is what
    ## `flat_tol = 1e-2` buys.  Two arms, and they must disagree --
    ## `flat_tol` moves the limit, the LADDER DEPTH does not.  docs-46
    ## predicted the opposite (one decade of limit per decade of ladder) and
    ## offered the falsification; this is it.
    def _limit(**kw):
        last = None
        for e in range(2, 22):
            g22 = 2.0 * 10.0 ** (-e)
            sig, inf = algebraic_conditioning(_e5_circuit(g22), **kw)
            if inf['verdict'] != 'well-conditioned':
                break
            if abs(sig - g22) > 0.1 * g22:
                break
            last = g22
        return last

    assert _limit(flat_tol=1e-2, decades=8) == 2e-7
    for dec in (12, 20, 30):
        assert _limit(flat_tol=1e-2, decades=dec) == 2e-7, dec

    ## ⚠⚠ GATE THE MECHANISM, NOT ITS CONSEQUENCE.  The arm above passes only
    ## because the FLOOR GUARD eats the extra probes -- so if anyone ever
    ## relaxes `floor_k`, the ladder-depth arm starts moving and this gate
    ## would fire on a change that IMPROVED the routine.  The actual invariant
    ## is the usable point count.  (docs-46 read this code and made exactly
    ## that objection to the first version of this gate.)
    counts = set()
    for dec in (8, 12, 20, 30):
        _s, inf = algebraic_conditioning(_e5_circuit(2e-7), decades=dec)
        counts.add(sum(inf['usable']))
    assert len(counts) == 1, counts

    ## ⚠ AND BOTH KNOBS BIND, ABOUT A DECADE EACH -- a single-cause story was
    ## wrong twice, once in each direction.  `flat_tol` fires; `floor_k`
    ## decides how many points the flatness test ever sees.
    assert _limit(flat_tol=1e0, decades=30) == 2e-8
    assert _limit(floor_k=1.0, decades=30) == 2e-8
    assert _limit(flat_tol=1e0, floor_k=1.0, decades=30) == 2e-9


def test_lte_grid_folds_a_driven_circuit_on_the_drives_own_period_boundaries():
    """The fold on a DRIVEN circuit (2026-09-21, found the day the fold went
    in): its sources are functions of absolute time, its edges sit at fixed
    phases of the drive, so the fold's boundaries are `k T` and phase 0 is
    the drive's t = 0 -- not a crossing of the fastest state, the oscillator
    rule.  Folded on a crossing, the switched-capacitor sampler's finest
    steps landed at phases 0.26 .. 0.31 of the drive with the switch edges
    (the steepest output) at 0.75; `fold=False` had hidden it because the
    default `tstab` = 200 T starts the window at the drive's phase 0.  Pinned:
    the ten finest fractions of the folded grid start within 0.05 T of the
    ten steepest phases of the converged output.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    fclk, cval, kb, temp = 100e3, 100e-12, 1.38e-23, 300.0
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
    cir = build()
    p = PSS(cir, method='radau', reltol=1e-6)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        fr, seed = p.lte_grid(T, x0=np.zeros(cir.n), reltol=1e-5)
    fr = np.asarray(fr, float)
    assert abs(p.lte_period / T - 1.0) < 1e-12          # the drive's, exactly
    q = PSS(build(), method='radau', reltol=1e-6)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        q.solve(period=T, timestep=T / len(fr), x0=seed, grid=fr, maxiterations=40)
    assert q.converged
    ts = np.asarray(q.waveform[0], float)
    X = np.delete(np.asarray(q.waveform[1], float), q.irefnode, axis=0)
    out = [str(n_) for n_ in cir.nodes].index('out')
    out = out if out < q.irefnode else out - 1
    ## the fold's OWN fractions (the solve's grid adds `_period_grid`'s
    ## opener ramp and closing step at the seam, phase 0 -- structural,
    ## the same on every folded grid, and not where the fold put its steps)
    starts = np.concatenate([[0.0], np.cumsum(fr)[:-1]])
    ph_fine = np.sort(starts[np.argsort(fr)[:10]] % 1.0)
    dv = np.abs(np.gradient(X[out], ts))
    ph_steep = np.sort(ts[np.argsort(dv)[-10:]] / T % 1.0)
    d = np.abs(ph_fine[:, None] - ph_steep[None, :])
    d = np.minimum(d, 1.0 - d)
    assert np.max(np.min(d, axis=1)) < 0.05, (ph_fine, ph_steep)   # was 0.45 off


def test_the_period_quadrature_breaks_its_spline_at_landed_events_and_reaches_fourth_order_there():
    """Item 2 of the non-uniform-grid list (2026-09-21): under landed events
    every period integral kept the TRAPEZOID -- a spline C^2 across a kink
    rings -- and that capped every method's harmonics and noise at second
    order on exactly the grids a clocked circuit gets.  `periodic_spline_weights`
    now takes `breaks`: a not-a-knot cubic spline PER SEGMENT between the
    event nodes (and node 0), never crossing a kink.

    (a) Exact ladder: `exp(sin 2 pi t) + 3 tri(t - te)` (integral I_0(1) +
    3/4), both corners landed on a 1 + 0.15 sin grid.  Pinned: the
    piecewise rule at least 30x below the trapezoid at N = 202 and falling
    at least 8x per doubling to 802 (measured -6.3e-10 / -2.1e-11 /
    -5.9e-13 against the trapezoid's -3.6e-7 / -9.3e-8 / -2.5e-8); on the
    smooth integrand alone, breaking only at node 0 costs a constant, not
    the order (-8.2e-10 vs +6.4e-11 at 200).
    (b) The circuit: radau on its own fold of the pulse-clocked sampler,
    ramp ends landed, the output's harmonics from `_period_quadrature`
    against a radau uniform-3200 reference: H1 5.2e-2 -> 5.3e-4 and H2
    3.1e-1 -> 2.1e-2 at 90 points (measured; H1 reaches the reference's
    1.3e-5 floor by 177).  Pinned: |H1| within 2e-3 and |H2| within 5e-2
    of the reference, and the weights differ from the trapezoid's.
    """
    import warnings as _w
    from scipy.special import i0
    from pycircuit.circuit.shooting import periodic_spline_weights

    ## (a)
    te = 0.3137
    exact = float(i0(1.0)) + 0.75
    err_t, err_p = {}, {}
    for N in (202, 402, 802):
        u = np.linspace(0.0, 1.0, N - 2, endpoint=False)
        t = u + 0.15 * np.sin(2 * np.pi * u) / (2 * np.pi)
        t = np.sort(np.concatenate([t, [te, (te + 0.5) % 1.0]]))
        x = (t - te) % 1.0
        y = np.exp(np.sin(2 * np.pi * t)) + 3.0 * np.where(x < 0.5, x, 1.0 - x)
        br = [int(np.argmin(np.abs(t - te))),
              int(np.argmin(np.abs(t - ((te + 0.5) % 1.0))))]
        h = np.diff(np.r_[t, t[0] + 1.0])
        err_t[N] = abs(float((0.5 * (h + np.roll(h, 1))) @ y) - exact)
        err_p[N] = abs(float(periodic_spline_weights(t, 1.0, br) @ y) - exact)
    assert err_t[202] / err_p[202] > 30.0, (err_t, err_p)
    assert err_p[202] / err_p[402] > 8.0 and err_p[402] / err_p[802] > 8.0, err_p
    u = np.linspace(0.0, 1.0, 200, endpoint=False)
    t = u + 0.15 * np.sin(2 * np.pi * u) / (2 * np.pi)
    y = np.exp(np.sin(2 * np.pi * t))
    e_seam = abs(float(periodic_spline_weights(t, 1.0, [0]) @ y) - float(i0(1.0)))
    assert e_seam < 1e-8, e_seam

    ## (b)
    circuit.default_toolkit = circuit.numeric
    T = 1e-5

    def solve(npts, grid=None, seed=None, reltol=1e-8):
        cir = _pulse_clocked_sampler(T)
        p = PSS(cir, method='radau', reltol=reltol)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / npts, grid=grid,
                    x0=np.zeros(cir.n - 1) if seed is None else seed,
                    maxiterations=100)
        assert p.converged
        io = [str(n_) for n_ in cir.nodes].index('out')
        io = io if io < p.irefnode else io - 1
        fp = p.factored_period()
        tm = np.asarray(fp.times, float)[:len(fp.steps)]
        X = np.asarray(p.waveform[1], float)
        X = X if X.shape[0] == cir.n - 1 else np.delete(X, p.irefnode, axis=0)
        v = X[io][:len(tm)]
        w = p._period_quadrature(fp)
        if w is None:                        # uniform: the index DFT
            w = np.full(len(tm), 1.0 / len(tm))
        H = np.array([np.sum(w * v * np.exp(-2j * np.pi * k * tm / T))
                      for k in range(3)])
        return cir, p, H, w, tm

    _c, _p, Href, _w0, _t0 = solve(3200, reltol=1e-10)
    cir = _pulse_clocked_sampler(T)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        fr, seed = PSS(cir, method='radau', reltol=1e-8).lte_grid(
            T, x0=np.zeros(cir.n), reltol=1e-5)
    _c, p, H, w, tm = solve(len(fr), grid=np.asarray(fr, float), seed=seed)
    assert p.event_times, 'the ramp ends must be landed for this to test anything'
    h = np.diff(np.r_[tm, tm[0] + T])
    assert np.max(np.abs(w - 0.5 * (h + np.roll(h, 1)) / T)) > 1e-3 * np.max(w)
    rel = np.abs(H - Href) / np.abs(Href)
    assert rel[1] < 2e-3 and rel[2] < 5e-2, rel


def test_the_recurrence_detector_picks_one_crossing_per_period_by_the_hint():
    """Item 2 of the non-uniform-grid list (2026-09-21): `_observed_period`
    took every rising crossing of the fastest-swinging state, so a state
    that crosses its midline more than once per period -- a strong third
    harmonic (three crossings), two pulses of different height (two) --
    spread the spacings past the 1 % rule and returned None: the fold
    was refused (warned, "no consistent recurrence") and the grid fell
    back to the single window.  `_crossing_chain` now walks back from the
    last crossing choosing, per hint period, the crossing nearest the
    expected time; a single-crossing state gives the same chain as
    before, bit for bit.  Pinned on synthetic waveforms with the period
    known: the two multi-crossing cases return the period to 1e-6 from
    hints 10 % either side, the single-crossing case is unchanged, and
    the fold's boundaries are one per period on the three-crossing state.
    """
    T = 2.0
    w = 2.0 * np.pi / T
    t = np.sort(np.concatenate([np.linspace(0.0, 40 * T, 8000),
                                np.linspace(0.0, 40 * T, 8000) + 0.0013]))
    ph = (t % T) / T
    cases = {
        'sin': np.sin(w * t),
        'sin + 1.2 sin 3': np.sin(w * t) + 1.2 * np.sin(3 * w * t),
        'two pulses': (0.2 * np.sin(w * t)
                       + 3.0 * np.exp(-((ph - 0.2) / 0.02) ** 2)
                       + 2.5 * np.exp(-((ph - 0.6) / 0.02) ** 2)),
    }
    for label, v in cases.items():
        xs = np.vstack([v, 0.1 * v])
        y = v - 0.5 * (v.max() + v.min())
        per_period = np.sum((y[:-1] < 0.0) & (y[1:] >= 0.0)) / 40.0
        for hint in (T, 1.1 * T, 0.9 * T):
            r = PSS._observed_period(t, xs, 5, hint)
            assert r is not None, (label, hint)
            assert abs(r / T - 1.0) < 1e-6, (label, hint, r)
        if label == 'sin':
            assert per_period == 1.0
        else:
            assert per_period >= 2.0, (label, per_period)
    ## the fold's boundaries come from the same chain
    v = cases['sin + 1.2 sin 3']
    y = v - 0.5 * (v.max() + v.min())
    tc = PSS._crossing_chain(t, y, T)
    d = np.diff(tc)
    assert len(tc) >= 30 and np.max(np.abs(d / T - 1.0)) < 1e-6, (len(tc), d[:5])


def test_the_folds_walk_ends_at_the_period_so_two_folds_of_one_orbit_agree_and_gear_no_longer_swings():
    """Item 3 of the non-uniform-grid list (2026-09-21): the "gear swing" on
    folded grids was the FOLD's alignment, not gear.  `_fold_periods` walked
    past the period by up to one coarse step and rescaled every fraction to
    sum 1, shifting every phase by up to that step times its phase: six
    folds of vdP mu = 10 (205 points, the same edge groups to the point)
    split into two clusters whose fine groups sat 0.007 / 0.020 T apart, and
    gear read +1409 / +1563 / +1489 ppm on one cluster against +273 .. +302
    on the other -- trbdf2 on the SAME grids +150 .. +173 against +44 .. +46,
    a 5x split for every second-order method.  The walk now ends exactly
    at the period (the remainder spread over the trailing steps within the
    controller's 2x growth): gear +26 / +18 / +24 / +16 / -44 / +30, trbdf2
    +12 / +10 / +13 / +10 / +4 / +12.  Pinned on two folds from different
    transient starts: their fine groups' phases agree to 8e-3 T (measured
    0.0005 / 0.0042; were 0.007 / 0.020 apart), no consecutive ratio exceeds 2 in either grid, and gear's
    period on each is within 120 ppm of the reference and the two within
    40 ppm of each other (were 1100 ppm apart).
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

    folds = {}
    for v0 in (2.0, 1.0):
        cir = vdp()
        p = PSS(cir, method='gear')
        xfull = np.zeros(cir.n)
        xfull[[str(n_) for n_ in cir.nodes].index('v')] = v0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            fr, seed = p.lte_grid(19.1, x0=xfull, reltol=1e-5)
        fr = np.asarray(fr, float)
        assert abs(fr.sum() - 1.0) < 1e-12
        r = fr[1:] / fr[:-1]
        assert np.max(np.maximum(r, 1.0 / r)) <= 2.0 + 1e-9, np.max(np.maximum(r, 1.0 / r))
        st = np.concatenate([[0.0], np.cumsum(fr)[:-1]])
        fine = np.flatnonzero(fr < 0.004)
        groups = np.split(fine, np.flatnonzero(np.diff(fine) > 3) + 1)
        q = PSS(vdp(), method='gear', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=p.lte_period, timestep=p.lte_period / len(fr), x0=seed,
                    maxiterations=80, break_events=False, grid=fr)
        assert q.converged
        folds[v0] = ([float(st[g[0]]) for g in groups], (q.period / T_REF - 1.0) * 1e6)
    (ga, ea), (gb, eb) = folds[2.0], folds[1.0]
    assert len(ga) == len(gb) == 2, (ga, gb)
    assert max(abs(x - y) for x, y in zip(ga, gb)) < 8e-3, (ga, gb)     # 0.0005 / 0.0042 measured
    assert abs(ea) < 120 and abs(eb) < 120 and abs(ea - eb) < 40, (ea, eb)


def test_state_events_become_newton_unknowns_and_land_the_grid_on_a_pwm_switching_instant():
    """The event half of B7 (2026-09-21/22, Andreas: "Do the events as a
    Newton unknown").  On the PWM loop every method was FIRST order on a
    uniform grid because the comparator's crossing sits inside a step
    (radau, mean error 2.9e-2 -> 3.3e-3 V over 100 -> 800 points, halving
    per doubling).  With `state_events` (the default) a driven radau or
    trbdf2 solve runs a bordered second stage: the crossings of the
    first stage's orbit (both edges of the switch's window, declared by
    `VSwitch.state_events()`) become unknowns with `x_0`, the grid between
    consecutive events scales with its segment, and the Newton lands the
    grid on the crossings.  ⚠ Two things the period column never needed
    were found by the FD check and are pinned here: a driven source moves
    with the grid (`-u_dot . dt_stage` on every stage), and `_k_at` takes
    the stage time (it took t = 0, "autonomous only").  Pinned: the
    bordered Jacobian's event column and row against central differences
    to 1e-6 (radau and trbdf2); the staged solve at 100 points at least
    5x closer to a staged 800-point radau reference than the
    unstaged one, with the on-crossing within 5e-4 of the
    reference's; `state_events=False` reproduces the one-stage solve;
    gear stages (phase B), and so do the plain one-step map opened at
    x(0) (trap) and a Nordsieck GLM (2026-09-24).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5

    def solve(method, N, se, reltol=1e-8):
        cir = _pwm_loop(T)
        p = PSS(cir, method=method, reltol=reltol)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            p.solve(period=T, timestep=T / N, x0=np.zeros(cir.n - 1),
                    maxiterations=100, state_events=se)
        assert p.converged
        io = [str(n_) for n_ in cir.nodes].index('out')
        io = io if io < p.irefnode else io - 1
        ts = np.asarray(p.waveform[0], float)
        X = np.asarray(p.waveform[1], float)
        v = X[io] if X.shape[0] == cir.n - 1 else np.delete(X, p.irefnode, axis=0)[io]
        return p, ts, v, [str(w_.message) for w_ in rec]

    ## (a) the Jacobian, against central differences, both kinds
    for method in ('radau', 'trbdf2'):
        p, ts, _v, _ws = solve(method, 60, False)
        x0 = np.asarray(p._period_state[1], float)
        m = p.cir.n - 1
        hs = np.diff(ts)
        W, c = p._state_event_rows()
        assert W is not None and W.shape[0] == 2
        base2, th0 = p._land_fractions(hs / T, [0.706, 0.7107])

        def FJ(z):
            xx, th = z[:m], z[m:]
            fr, hsens, nodes = p._event_remap(base2, th0, th, T)
            hs_ = fr * T
            tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
            wk = p._walk('stage', xx, tms_, hs_, T=T, hsens=hsens, capture=set(nodes))
            x_end, Mx, Pk = wk.x_end, wk.P, wk.Pk
            F = np.concatenate((xx - np.asarray(x_end),
                                [float(W[k] @ p._captured[nd][0]) - c[k]
                                 for k, nd in enumerate(nodes)]))
            J = np.zeros((m + 2, m + 2))
            J[:m, :m] = np.eye(m) - Mx
            for k in range(2):
                J[:m, m + k] = -np.asarray(Pk[k]).ravel()
            for k, nd in enumerate(nodes):
                xj, Pj, Pkj = p._captured[nd]
                J[m + k, :m] = W[k] @ Pj
                for l in range(2):
                    J[m + k, m + l] = float(W[k] @ Pkj[l])
            return F, J
        z = np.concatenate((x0, th0 + np.array([0.001, -0.001])))
        _F0, J0 = FJ(z)
        for i in (m, m + 1):
            zp = z.copy(); zp[i] += 1e-6
            zm = z.copy(); zm[i] -= 1e-6
            fd = (FJ(zp)[0] - FJ(zm)[0]) / 2e-6
            assert np.linalg.norm(J0[:, i] - fd) < 1e-6 * np.linalg.norm(fd), (method, i)

    ## (b) the accuracy, against a staged reference
    pr, tsr, vr, _ = solve('radau', 800, True, 1e-10)
    ref = lambda t: np.interp(t % T, tsr, vr)
    swing = vr.max() - vr.min()
    errs = {}
    for se in (False, True):
        p, ts, v, ws = solve('radau', 100, se)
        errs[se] = np.max(np.abs(v - ref(ts))) / swing
        if se:
            th = np.asarray(p._state_event_fracs, float)
            assert len(th) == 4, th
            assert abs(th[0] - pr._state_event_fracs[0]) < 5e-4, (th, pr._state_event_fracs)
            assert all(float(t_) in [float(e) for e in p.event_times] for t_ in th)
        else:
            assert p._state_event_fracs is None
    assert errs[False] / errs[True] > 5, errs

    ## (c) gear stages (phase B, 2026-09-22), the plain map opened at x(0)
    ## and a Nordsieck GLM too (2026-09-24)
    for meth in ('gear', 'trap', 'glm2'):
        _g, _ts, _v, _ws = solve(meth, 100, True)
        assert _g._state_event_fracs is not None and len(_g._state_event_fracs) == 4, meth


def test_the_staged_solves_monodromy_is_the_total_derivative_through_the_moving_event():
    """Phase B of events-as-unknowns (2026-09-22): the period map's
    derivative is NOT the fixed-grid monodromy.  A perturbation of `x_0`
    moves the crossing (`dtheta/dx_0 = -Gt^-1 G` from the event rows) and
    the state at the period moves with it through the event columns; the
    bordered system's Schur complement `M + P_theta dtheta/dx_0` is the
    saltation matrix, derived.  Measured on the PWM loop at 100 points:
    the fixed-grid `Mx` 109 % off the finite difference of the staged
    period map, the total 3.3e-8; dominant multiplier 0.691 where `Mx`
    reads 0.632 (and a LOCAL two-step saltation 0.656 -- which steps
    absorb the event's motion is not a small choice at this resolution,
    so the consumers need the bordered system, not a per-step patch).
    Pinned at 60 points: two columns of `_monodromy` against central
    differences of the staged map (inner Newton on theta) to 1e-5, and
    the fixed grid at least 30 % off.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    cir = _pwm_loop(T)
    p = PSS(cir, method='radau', reltol=1e-8)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1), maxiterations=100)
    assert p.converged and p._state_event_fracs is not None
    x0 = np.asarray(p._period_state[1], float)
    m = cir.n - 1
    base2 = np.asarray(p._grid_fracs, float)
    th_s = np.asarray(p._state_event_fracs, float)
    K = len(th_s)
    W, c = p._state_event_rows()
    Mtot = np.asarray(p._monodromy, float)
    assert np.asarray(p._event_sensitivity).shape == (K, m)

    def pieces(xx, th):
        fr, hsens, nodes = p._event_remap(base2, th_s, th, T)
        hs_ = fr * T
        tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
        wk = p._walk('stage', xx, tms_, hs_, T=T, hsens=hsens, capture=set(nodes))
        x_end, Mx = wk.x_end, wk.P
        gv = np.zeros(K)
        Gt = np.zeros((K, K))
        for k, nd in enumerate(nodes):
            xj, _Pj, Pkj = p._captured[nd]
            r = int(np.argmin([abs(float(W[i] @ xj) - c[i]) for i in range(W.shape[0])]))
            gv[k] = float(W[r] @ xj) - c[r]
            for l in range(K):
                Gt[k, l] = float(W[r] @ np.asarray(Pkj[l]).ravel())
        return np.asarray(x_end, float), np.asarray(Mx, float), gv, Gt

    def phi(xx):
        th = th_s.copy()
        for _it in range(30):
            x_end, Mx, gv, Gt = pieces(xx, th)
            if np.max(np.abs(gv)) < 1e-13:
                break
            th = th - np.linalg.solve(Gt, gv)
        return x_end, Mx

    _xe, Mx = phi(x0)
    ## the two columns that move the crossing most (the filtered output and the output)
    names = [str(n_) for i, n_ in enumerate(cir.nodes) if i != p.irefnode]
    for nm in ('fb', 'out'):
        i = names.index(nm)
        d = 1e-6 * max(1.0, abs(x0[i]))
        xp = x0.copy(); xp[i] += d
        xm = x0.copy(); xm[i] -= d
        col = (phi(xp)[0] - phi(xm)[0]) / (2 * d)
        assert np.linalg.norm(Mtot[:, i] - col) < 1e-5 * np.linalg.norm(col), nm
        ## with VSwitch's compact transition (2026-09-22) the fixed-grid map
        ## through the landed, resolved window is 8 % off the total derivative
        ## on `fb` (109 % with the tanh, whose tails sat outside the window);
        ## with the window in 16 steps (2026-09-25) 0.53 % -- part of the 8 %
        ## was the 8-step window's own sensitivity to where it sat
        assert np.linalg.norm(Mx[:, i] - col) > 0.003 * np.linalg.norm(col), nm



def test_a_staged_solve_reports_the_landed_map_on_every_kind(monkeypatch):
    """After a state-event stage, `_monodromy` -- and so `spectral_radius`
    -- is the map at the stage's final state on the LANDED grid: the total
    map through the events when their columns are built, the grid-frozen
    map when they cannot be (the one the unbordered consumers then use).

    ⚠ It is written once, after the stage.  Before, the stage kind wrote a
    partial map from inside the stage's Newton and gear's pair wrote none,
    so with the column assembly failing gear reported STAGE 1's map, from
    another grid and another orbit: measured on this loop, spectral radius
    0.905 against the landed map's 0.730.  (With the columns built, both
    kinds already reported the total map.)
    """
    import warnings as _w
    from pycircuit.circuit.shooting import events as _events
    circuit.default_toolkit = circuit.numeric
    T = 1e-5

    def solve(method):
        cir = _pwm_loop(T)
        p = PSS(cir, method=method, reltol=1e-8)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            p.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1),
                    maxiterations=100, state_events=True)
        assert p.converged, method
        fp = p.factored_period()
        Md = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                              for e in np.eye(fp.width)])
        return p, Md, [str(r.message) for r in rec]

    def rho(M):
        return float(np.max(np.abs(np.linalg.eigvals(np.asarray(M, dtype=float)))))

    for method in ('radau', 'gear'):
        p, Md, _msgs = solve(method)
        assert p._event_columns is not None, method
        Mt = Md + (np.asarray(p._event_columns['P_end'], dtype=float)
                   @ np.asarray(p._event_sensitivity, dtype=float))
        assert np.max(np.abs(np.asarray(p._monodromy) - Mt)) < 1e-8 * np.max(np.abs(Mt)), method
        assert abs(p.spectral_radius / rho(Mt) - 1.0) < 1e-9, (method, p.spectral_radius)

    def refuse(*a, **k):
        raise ValueError('forced')
    monkeypatch.setattr(_events.EventColumns, 'from_capture', staticmethod(refuse))
    for method in ('radau', 'gear'):
        p, Md, msgs = solve(method)
        assert p._event_columns is None, method
        assert any('could not assemble its event columns' in m for m in msgs), (method, msgs)
        assert np.max(np.abs(np.asarray(p._monodromy) - Md)) < 1e-8 * np.max(np.abs(Md)), method
        assert abs(p.spectral_radius / rho(Md) - 1.0) < 1e-9, (method, p.spectral_radius, rho(Md))


def test_an_autonomous_solve_takes_its_state_events_and_its_period_as_unknowns_together():
    """Phase B of events-as-unknowns (2026-09-22): on an AUTONOMOUS circuit
    the period joins the event unknowns, `z = [x_0, theta, T]`, the period
    one more column of the event algebra (`d h_j / d T = fraction_j`), the
    phase row closing the system.  The comparator relaxation oscillator
    with a sharp window: unstaged radau on uniform grids reads a period
    error of +4.4e-3 / +2.1e-4 / -1.7e-3 / -7.7e-4 at 100 / 200 / 400 / 800
    points against an unstaged 3200-point solve (itself untrusted below
    1e-4) -- set by where the crossing falls in its step, not by N.
    Against a STAGED radau-1600 reference the staged period reads
    -6.4e-4 / -1.9e-4 / -6.8e-5 / +8.1e-7 at 100 / 200 / 400 / 800 points
    (a thousandfold at 800, all four edges at the reference's).  ⚠ FD
    lesson: a theta step of 1e-6 across a 0.04 ns window is 3.5 % of it
    and reads the first-edge columns 1e-4 off by truncation; at 1e-8 every
    column is below 1e-7.  Pinned (measured -6.4e-4 vs -1.9e-4 staged,
    +4.4e-3 vs +2.0e-4 unstaged): the staged periods at 100 and 200
    points within 1e-3 of each other, the unstaged ones more than 2e-3
    apart and the unstaged more than 2e-3 from the staged at 100; four
    event fractions landed.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir = _comparator_relaxation_oscillator()
    seed, Tl = _relaxation_oscillator_seed(cir)

    def solve(N, se):
        c2 = _comparator_relaxation_oscillator()
        q = PSS(c2, method='radau', reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=Tl, timestep=Tl / N, x0=seed, maxiterations=100, state_events=se)
        assert q.converged
        return q

    q1, q2 = solve(100, True), solve(200, True)
    u1, u2 = solve(100, False), solve(200, False)
    assert q1._state_event_fracs is not None and len(q1._state_event_fracs) >= 2
    assert abs(q1.period / q2.period - 1.0) < 1e-3, (q1.period, q2.period)
    assert abs(u1.period / u2.period - 1.0) > 2e-3, (u1.period, u2.period)
    assert abs(u1.period / q1.period - 1.0) > 2e-3
    assert len(q1._state_event_fracs) == 4



def test_gear_lands_the_state_events_of_an_autonomous_orbit():
    """Gear on the comparator relaxation oscillator (2026-09-24, Andreas:
    "Analyse, fix and create a test").  Its free-period pair was the one kind
    the state-event stage skipped -- silently, since the method check counts
    gear as event-capable -- so it returned the UNSTAGED orbit: the crossing
    inside a step, the period set by where it falls there (+5.8e-3 of the
    exact period at 200 points, -2.9e-3 at 100, +3.6e-3 at 300), and the
    fixed-grid map's dominant multiplier 52.8 in place of 1.  Staged:
    +1.1e-4 at 200 points (+8.7e-5 / -6.8e-4 / -1.9e-4 / +1.3e-5 at 150 /
    100 / 300 / 800), a unit multiplier, the four window edges landed.  What
    is left is gear's own second order on the 10 ns discharge after the
    switch closes (its crossings trail radau's by 5e-3 T at 200 points,
    halving per doubling); radau lands 9e-8.

    Two false warnings the staged grid then raised, pinned silent at 100
    points.  The MULTIPLE-of-the-fundamental detector excluded its edges in
    POINTS, and gear's grid opens with ten doubling steps from 1e-5 T, so
    the orbit still leaving `x_0` read as a 3182-fold recurrence.  And the
    step-ratio warning counted a landed window's exit -- its ramp, then the
    partial step back to the base grid -- as REPEATED up-steps; smoothing
    those pairs away never improved the period (100, 150, 200, 800 points).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    ## `Tl` is the exact model's period
    seed, Tl = _relaxation_oscillator_seed(_comparator_relaxation_oscillator())

    def solve(N, se):
        q = PSS(_comparator_relaxation_oscillator(), method='gear', reltol=1e-9)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=Tl, timestep=Tl / N, x0=seed, maxiterations=100,
                    state_events=se)
        assert q.converged, (N, se)
        return q, [str(r.message) for r in rec]

    def err(q):
        return abs(q.period / Tl - 1.0)

    def rho(q):
        return float(np.max(np.abs(q.floquet_multipliers)))

    staged, msgs = solve(200, True)
    unstaged, _m = solve(200, False)
    assert len(staged._state_event_fracs) == 4 and staged.stage_one_converged is True
    ## ⚠ RE-PINNED 2026-09-25 for the 16-step window (`event_window_steps`):
    ## 3.0e-4 at 200 points.  The 1.1e-4 at 8 steps was a CANCELLATION -- the
    ## 8-step ladder alternated in sign (-6.3e-4 / +1.1e-4 / -1.2e-4 / +1.3e-5
    ## at 100 / 200 / 400 / 800); at 16, -2.2e-4 / -3.0e-4 / -1.2e-4 / +1.4e-5,
    ## the two agreeing from 400 on, where the window no longer dominates
    assert err(staged) < 5e-4 and err(unstaged) > 3e-3, (err(staged), err(unstaged))
    assert abs(rho(staged) - 1.0) < 1e-3 and rho(unstaged) > 10.0, (rho(staged), rho(unstaged))
    q100, msgs100 = solve(100, True)
    assert len(q100._state_event_fracs) == 4 and err(q100) < 1.5e-3, err(q100)
    for m in msgs + msgs100:
        assert 'MULTIPLE' not in m and 'steps up by' not in m, m


def test_a_state_event_stage_rescues_a_first_stage_that_did_not_converge():
    """Across a sharp switch the UNSTAGED map is nearly discontinuous in the
    state -- where the crossing falls in its step moves with every iterate
    -- and its Newton can fail where the staged system converges.  Measured
    on the comparator oscillator under gear: the unstaged solve fails at 350
    and 400 points (the fixed-grid map's dominant multiplier 100-290, and a
    'solved-history stall' diagnosed although the second multiplier is
    0.002), and the stage, run from its last iterate, converges to 1.5e-4 /
    1.2e-4 of the exact period.

    The stage's own Newton decides `converged`; `stage_one_converged`
    records that the first stage did not; and the first stage's stall
    diagnosis -- about a system the result no longer comes from -- is held
    back.  When the stage fails too (15 iterations leave it short) that
    diagnosis is emitted after all, and it names the state events.  30
    iterations keep the rescue quick.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    seed, Tl = _relaxation_oscillator_seed(_comparator_relaxation_oscillator())

    def solve(maxiterations):
        q = PSS(_comparator_relaxation_oscillator(), method='gear', reltol=1e-9)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=Tl, timestep=Tl / 350, x0=seed,
                    maxiterations=maxiterations, state_events=True)
        return q, [str(r.message) for r in rec]

    q, msgs = solve(30)
    assert q.converged and q.stage_one_converged is False
    assert len(q._state_event_fracs) == 4
    assert abs(q.period / Tl - 1.0) < 3e-4, q.period / Tl - 1.0
    assert not any('stall' in m or 'did not converge' in m for m in msgs), msgs

    q, msgs = solve(15)
    assert not q.converged and q.stage_one_converged is False
    stall = [m for m in msgs if 'solved-history stall' in m]
    assert len(stall) == 1 and 'state_events=True' in stall[0], msgs



def test_the_plain_map_lands_state_events_opened_at_x0():
    """The state-event stage skipped the PLAIN map (euler, trap, theta):
    `_walk_lmm` refused its event columns.  Opened as the plain map is by
    default -- `x(0)` manufactured by an order-dropped step -- the stage does
    not work: that step moves with the first crossing, and on trap it failed
    to converge (comparator oscillator, 200 points) or read worse than the
    unstaged solve (800).  OPENED AT `x(0)` (`x0_unknown`, now its default
    when the circuit declares events) it does -- after two defects in the
    columns, both found against finite differences on a driven PWM loop:

    * the source's motion (`u_dot`) went into the CARRIED companion
      sensitivity `Pq`; it is a current, not a charge.  Gear's pair never
      reads `Pq` (`b = 0`), trap's does -- its columns were 0.3-380x off;
    * the previous step's partial, estimated from the total under uniform
      scaling, assumes coefficients homogeneous in `h`; theta's are not, and
      on a one-step companion that partial is zero by structure (theta's
      columns up to 358x off).

    Measured after: trap on the oscillator +1.4e-4 at 200 points against
    +6.4e-3 unstaged; on the PWM loop, where the unstaged trap switches at
    0.51 of the period against 0.688, the staged crossings within 1.2e-4 of
    radau's (400 points).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    seed, Tl = _relaxation_oscillator_seed(_comparator_relaxation_oscillator())

    def osc(method, se):
        q = PSS(_comparator_relaxation_oscillator(), method=method, reltol=1e-9)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=60,
                    state_events=se)
        assert q.converged, (method, se)
        return q, [str(r.message) for r in rec]

    staged, msgs = osc('trap', True)
    unstaged, _m = osc('trap', False)
    assert staged._open_at_x0 and len(staged._state_event_fracs) == 4
    assert abs(staged.period / Tl - 1.0) < 5e-4 < 3e-3 < abs(unstaged.period / Tl - 1.0), \
        (staged.period / Tl - 1.0, unstaged.period / Tl - 1.0)
    assert not any('state event' in m for m in msgs), msgs

    ## the driven loop: the crossings against radau's, and the columns
    T = 1e-5
    th_radau = np.array([0.688459, 0.693154, 0.992922, 0.992972])   # radau, 400 pts
    p = PSS(_pwm_loop(T), method='trap', reltol=1e-10)
    ## ⚠ SETTLED FIRST (`tstab`) since a landed ramp edge drops the order
    ## (2026-09-28): from zeros the staged Newton then wandered (337
    ## evaluations, no convergence), though from a settled seed it converges
    ## in 14 to the same crossings, its event columns FD-exact -- a basin,
    ## not a Jacobian (the UNSTAGED solve converged from zeros only WITH the
    ## drop).  `_staged_fallback` recovers the zero start since, at ~9x the
    ## time (`test_a_staged_solve_that_stalls_falls_back_to_the_one_stage_orbit`).
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T, timestep=T / 60, x0=np.zeros(_pwm_loop(T).n - 1),
                maxiterations=100, tstab=20 * T)
    assert p.converged and p._open_at_x0
    assert np.max(np.abs(np.asarray(p._state_event_fracs) - th_radau)) < 1e-3

    for method in ('trap', 'theta'):
        q = PSS(_pwm_loop(T), method=method, reltol=1e-10)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=T, timestep=T / 60, x0=np.zeros(_pwm_loop(T).n - 1),
                    maxiterations=100, state_events=False)
        x0 = np.asarray(q._period_state[1], dtype=float)[:q.cir.n - 1]
        W, c = q._state_event_rows()
        times, hs = q._period_grid(T, 60, None)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q._walk('plain', x0, times, hs, T=T, hsens=np.zeros((len(hs), 0)),
                    capture=set(range(1, len(hs) + 1)), open_at_x0=True)
        base2, th0, Wk, ck = q._stage_one_crossings(x0, times, T, W, c)

        def end_and_cols(th):
            fr, hsens, nodes = q._event_remap(base2, th0, th, T)
            hs_ = fr * T
            tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                w = q._walk('plain', x0, tms_, hs_, T=T, hsens=hsens,
                            capture=set(nodes), open_at_x0=True)
            return (np.asarray(w.x_end, dtype=float),
                    np.column_stack([pk[0] for pk in w.Pk]))
        _xe, Pk = end_and_cols(th0)
        ## ⚠ 1e-7, not 1e-8: with the 16-step window (sub-steps 5e-6 T) the
        ## difference is roundoff-limited below it -- theta's third column
        ## read 1.6e-4 / 1.8e-3 / 1.6e-1 at eps 1e-7 / 1e-8 / 1e-9, growing
        ## as eps shrinks (measured 2026-09-25)
        eps = 1e-7
        for l in range(len(th0)):
            d = np.zeros(len(th0))
            d[l] = eps
            fd = (end_and_cols(th0 + d)[0] - end_and_cols(th0 - d)[0]) / (2 * eps)
            err = np.max(np.abs(Pk[:, l] - fd)) / np.max(np.abs(fd))
            assert err < 1e-3, (method, l, err)

    ## an explicit x0_unknown=False keeps the opener, and the stage says so
    q = PSS(_comparator_relaxation_oscillator(), method='trap', reltol=1e-9)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        q.solve(period=Tl, timestep=Tl / 200, x0=seed, maxiterations=60,
                x0_unknown=False)
    assert q._state_event_fracs is None
    assert any('needs the map opened at x(0)' in str(r.message) for r in rec)


def test_a_state_event_stage_that_fails_hands_back_the_first_stage():
    """A stage iterate can leave the orbit far enough that an inner step
    does not converge -- measured: trap's first event columns on a driven
    PWM loop, before they were right, sent `vin` to 1e9 -- and the whole
    solve RAISED.  The stage improves on stage 1; it is not a condition of
    having a result: a failure warns and returns the first stage.
    """
    import types
    import warnings as _w
    from pycircuit.circuit.analysis import NoConvergenceError
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    p = PSS(_pwm_loop(T), method='radau', reltol=1e-9)

    def boom(self, *a, **k):
        raise NoConvergenceError('forced')
    p._event_remap = types.MethodType(boom, p)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        p.solve(period=T, timestep=T / 60, x0=np.zeros(_pwm_loop(T).n - 1),
                maxiterations=100)
    assert p.converged and p._state_event_fracs is None
    assert any('state-event stage failed' in str(r.message) for r in rec)

def test_gears_state_event_stage_carries_both_step_partials_and_lands_the_crossing():
    """Phase B of events-as-unknowns (2026-09-22): gear's event columns.  A
    two-step companion's residual depends on its own step AND the previous
    one, so each column carries `dr/dh_n hsens[j] + dr/dh_{n-1} hsens[j-1]`
    with the previous step's partial assembled from the two that exist --
    `residual_dT` is `sum_j h_j dr/dh_j` (Euler's theorem, no source) and
    `residual_dh` is the coefficients' partial PLUS `du/dt`, so `dr/dh_{n-1}
    = (residual_dT - (residual_dh - u_dot) h_n) / h_{n-1}` -- and the
    source's motion once: `residual_dh`'s own part and `u_dot tau_n` for
    the shift of the step's start.  ⚠ Counting `u_dot (tau + w)` on top
    of `residual_dh` read the node after a landed event 153 % off.  Pinned:
    the columns against central differences to 1e-6 at the period and at
    the event nodes (PWM loop, two anchors), and the staged gear solve at
    200 points converging with its crossing within 5e-4 T of the staged
    radau reference's and no worse than the unstaged one (measured 4.0e-3
    vs 4.4e-3 of the swing; gear's own second-order off-phase error
    dominates here, as trbdf2's does).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5

    def solve(method, N, se):
        cir = _pwm_loop(T)
        p = PSS(cir, method=method, reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / N, x0=np.zeros(cir.n - 1),
                    maxiterations=100, state_events=se)
        assert p.converged
        io = [str(n_) for n_ in cir.nodes].index('out')
        io = io if io < p.irefnode else io - 1
        ts = np.asarray(p.waveform[0], float)
        X = np.asarray(p.waveform[1], float)
        v = X[io] if X.shape[0] == cir.n - 1 else np.delete(X, p.irefnode, axis=0)[io]
        return p, ts, v

    ## (a) the columns
    p, ts, _v = solve('gear', 100, False)
    x0s = np.asarray(p._period_state[1], float)
    m = p.cir.n - 1
    hs = np.diff(ts)
    X = np.delete(np.asarray(p.waveform[1], float), p.irefnode, axis=0)
    xm1 = X[:, -2].copy()
    base2, th0 = p._land_fractions(hs / T, [0.6914, 0.6961])
    K = 2

    def run(th):
        fr, hsens, nodes = p._event_remap(base2, th0, th, T)
        hs_ = fr * T
        tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
        wk = p._walk('pair', np.concatenate((x0s, xm1)), tms_, hs_, T=T,
                     hsens=hsens, capture=set(nodes))
        return (np.asarray(wk.x_end), [np.asarray(pk[0]).ravel() for pk in wk.Pk], nodes,
                {nd: (np.asarray(p._captured[nd][0]), [np.asarray(c).ravel() for c in p._captured[nd][2]]) for nd in nodes})

    th = th0 + np.array([0.002, -0.003])
    xl, Pkl, nodes, caps = run(th)
    for k in range(K):
        d = 1e-6
        tp = th.copy(); tp[k] += d
        tm = th.copy(); tm[k] -= d
        xlp, _, _, cp = run(tp)
        xlm, _, _, cm = run(tm)
        fd = (xlp - xlm) / (2 * d)
        assert np.linalg.norm(Pkl[k] - fd) < 1e-6 * np.linalg.norm(fd), k
        for nd in nodes:
            fdn = (cp[nd][0] - cm[nd][0]) / (2 * d)
            assert np.linalg.norm(caps[nd][1][k] - fdn) < 5e-6 * np.linalg.norm(fdn), (k, nd)   # 1.4e-6 measured at one node: the FD's own noise

    ## (b) the stage
    pr, tsr, vr = solve('radau', 800, True)
    ref = lambda t: np.interp(t % T, tsr, vr)
    swing = vr.max() - vr.min()
    qs, tss, vs = solve('gear', 200, True)
    qu, tsu, vu = solve('gear', 200, False)
    assert qs._state_event_fracs is not None and len(qs._state_event_fracs) == 4
    assert abs(qs._state_event_fracs[0] - pr._state_event_fracs[0]) < 5e-4
    es = np.max(np.abs(vs - ref(tss))) / swing
    eu = np.max(np.abs(vu - ref(tsu))) / swing
    ## ⚠ with VSwitch's compact transition (2026-09-22) gear's stage buys
    ## NOTHING on this loop: staged 4.3e-3 against unstaged 3.7e-3 at 200
    ## points -- like trbdf2's, gear's own second-order off-phase error is
    ## what is left once the transition sits inside the landed window.  The
    ## columns are exact (above); the pin is that the stage does no harm
    ## beyond that margin.  Radau is the method for a staged solve.
    assert es <= 1.3 * eu, (es, eu)


def _windowed_exact_relaxation_period():
    """The EXACT period of `_comparator_relaxation_oscillator` WITH its
    0.2 mV window: `VSwitch`'s compact smoothstep conductance integrated at
    rtol 1e-12 (DOP853), the period by shooting on the falling crossing of
    `fb1 - 2.5`.  -7.997e-8 off the ideal-comparator period; the reference
    for the staged solve's fifth order AT the switch."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve
    R1 = R2 = R3 = 1e3
    C1, C2, C3 = 1e-9, 3e-10, 3e-10
    VDD, VREF, RON, ROFF, VON, VOFF = 5.0, 2.5, 10.0, 1e7, 1e-4, -1e-4

    def g_of(vc):
        t = min(max((vc - VOFF) / (VON - VOFF), 0.0), 1.0)
        return 1 / ROFF + (1 / RON - 1 / ROFF) * t * t * t * (10 + t * (-15 + 6 * t))

    def f(t, x):
        c, fb0, fb1 = x
        g = g_of(fb1 - VREF)
        return [((VDD - c) / R1 - (c - fb0) / R2 - g * c) / C1,
                ((c - fb0) / R2 - (fb0 - fb1) / R3) / C2,
                ((fb0 - fb1) / R3) / C3]

    def cross(t, x):
        return x[2] - VREF
    cross.direction = -1
    cross.terminal = True

    def one_period(x0):
        a = solve_ivp(f, (0, 2e-7), x0, method='DOP853', rtol=1e-12, atol=1e-14, max_step=2e-9)
        s = solve_ivp(f, (2e-7, 5e-6), a.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-14,
                      events=cross, max_step=2e-9)
        assert len(s.t_events[0]) == 1
        return s.y_events[0][0], float(s.t_events[0][0])

    def resid(z):
        x1, _T = one_period([z[0], z[1], VREF])
        return [x1[0] - z[0], x1[1] - z[1]]

    z = fsolve(resid, [0.0717, 2.2036], xtol=1e-12)
    _x1, T = one_period([z[0], z[1], VREF])
    return T


def test_the_staged_solve_is_fifth_order_at_the_switch_against_the_windowed_exact_period():
    """E8 closed (2026-09-22): with `VSwitch`'s compact transition the
    staged radau period converges to the WINDOWED exact period as -7.3e-7 /
    -1.85e-8 / -9.0e-10 / -2.5e-10 at 100 / 200 / 400 / 800 points (with
    the tanh it was -6.7e-4 / -2.2e-4 / -9.6e-5 / -2.7e-5 against the ideal
    one and never better than first order: 12 % of the transition sat
    outside the landed window).  Pinned at 100 and 200 points: within 2e-6
    and 5e-8 of the windowed exact period, and the ratio between them
    above 8 (fifth order would be 32; the 100-point solve is where the
    grid's own resolution of the two lags still shows)."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T_ex = _windowed_exact_relaxation_period()
    assert abs(T_ex / 1.3918374e-6 - 1.0) < 2e-7, T_ex
    cir = _comparator_relaxation_oscillator()
    names = [str(n_) for n_ in cir.nodes]
    seed, Tl = _relaxation_oscillator_seed(cir)
    errs = {}
    ## ⚠ RE-PINNED 2026-09-25 for the 16-step window (`event_window_steps`):
    ## against the windowed exact period, -5.6e-9 / -6.0e-9 / -6.3e-10 /
    ## -2.7e-11 at 100 / 200 / 400 / 800 (at 8 steps -1.8e-7 / -1.8e-8 /
    ## -8.5e-10 / -2.5e-10) -- lower everywhere, and the order shows from
    ## 200 on (9.5, 23): at 100 points the window no longer dominates
    for N in (200, 400):
        q = PSS(cir, method='radau', reltol=1e-10)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=Tl, timestep=Tl / N, x0=seed, maxiterations=100, state_events=True)
        assert q.converged
        errs[N] = abs(q.period / T_ex - 1.0)
    assert errs[200] < 2e-8 and errs[400] < 2e-9, errs
    assert errs[200] / errs[400] > 8.0, errs


def test_gears_stage_stores_its_event_columns_and_its_bordered_pac_matches_radaus():
    """Item 3 of the 2026-09-22 list: the gear stage stored no event columns,
    so every bordered consumer ran unbordered on a staged gear solve.  Now it
    stores them in gear's PAIR form (`P_nodes[j]` m x 2m, `P_end` 2m x K,
    the total monodromy the pair map plus the saltation) and `PAC.solve`'s
    bordered block takes the pair width.  Measured on the PWM loop at 100
    points against radau's bordered response (exact against the fixed-time
    FD): gear bordered 5.6e-4 / 2.5e-4 off (out / fb), gear UNBORDERED
    10.6 % -- a two-step map through the landed window does not carry the
    switching the way radau's one-step map does, so here the bordering is
    worth two orders of magnitude, not a sliver.  Gear's pair-total
    multipliers 0.98106 / 0.68206 against radau's 0.98105 / 0.68203.  The
    covariance closure and the adjoint row stay unbordered on gear (warned,
    plain)."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6

    def staged(method):
        cir = _pwm_loop(T)
        p = PSS(cir, method=method, reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            p.solve(period=T, timestep=T / 100, maxiterations=100, state_events=True)
        assert p.converged
        return cir, p

    cr, pr = staged('radau')
    cg, pg = staged('gear')
    assert pg._event_columns is not None and np.shape(pg._monodromy) == (2 * (cg.n - 1), 2 * (cg.n - 1))
    lr = np.sort(np.abs(np.linalg.eigvals(np.asarray(pr._monodromy))))[::-1]
    lg = np.sort(np.abs(np.linalg.eigvals(np.asarray(pg._monodromy))))[::-1]
    assert np.max(np.abs(lg[:2] - lr[:2]) / lr[:2]) < 1e-3, (lg[:2], lr[:2])
    names = [str(n_) for n_ in cr.nodes]
    red = [nm for i, nm in enumerate(names) if i != pr.irefnode]
    f0 = 1.0 / T

    def resp(cir, p, bordered):
        ev = p._event_columns
        if not bordered:
            p._event_columns = None
        try:
            pac = PAC(cir, toolkit=circuit.numeric)
            with _w.catch_warnings():
                _w.simplefilter('ignore')
                pac.solve(p, [f0])
            tt, yy = pac.time_response[0]
        finally:
            p._event_columns = ev
        return np.asarray(tt, dtype=float), np.asarray(yy)

    tr_, yr = resp(cr, pr, True)
    for bordered, lo, hi in ((True, 0.0, 2e-3), (False, 5e-2, 1.0)):
        tg_, yg = resp(cg, pg, bordered)
        for nm in ('out', 'fb'):
            i = red.index(nm)
            yref = (np.interp(tg_, tr_, np.real(yr[:, i]))
                    + 1j * np.interp(tg_, tr_, np.imag(yr[:, i])))
            err = float(np.max(np.abs(yg[:, i] - yref)) / np.max(np.abs(yref)))
            assert lo < err < hi, (bordered, nm, err)


def test_a_glm_lands_state_events_and_restarts_where_its_step_grows():
    """The state-event stage skipped a Nordsieck GLM (a warning: the
    crossing inside a step, first order there).  Built: the GLM walk
    carries event columns -- each step's own `h` in its stage and output
    rows, the entering rescale's ``d rho``, a driven source moving with the
    stage times, the startup's own substeps at node 0 -- and the stage
    reads the map on the state as it does a stage method's.

    ⚠ THE LANDED GRID NEEDED A RESTART.  Its window between a switch's two
    edges is a few short steps, and leaving it the step grows up to 1125x:
    rescaling the Nordsieck vector by `rho^k` extrapolates the transition's
    derivatives over the long step.  On that grid alone (no stage) glm3
    read 1.2e-1 of the PWM loop's swing off radau and glm2 1.5e-3, and the
    stage converged to a false crossing (0.8685 against 0.688) or failed.
    The transient now restarts where the step grows more than
    `GLM_RESTART_GROWTH` (4x; its adaptive controller grows at most 2x),
    and the linearisation follows: the step before ends on the startup of
    its last stage.

    Pinned: the bordered Jacobian against central differences on the PWM
    loop (glm2, with a restart in the walk); on the comparator relaxation
    oscillator against its exact period, glm3 staged +2.5e-7 at 100 points
    where unstaged reads -3.6e-3 (and a fixed-grid multiplier of 0.30 for
    the unit one), glm2 staged +6.4e-5 at 200 (unstaged +6.0e-4).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    cir = _pwm_loop(T)
    p = PSS(cir, method='glm2', reltol=1e-8)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        p.solve(period=T, timestep=T / 60, x0=np.zeros(cir.n - 1),
                maxiterations=100, state_events=False)
    assert p.converged
    ts = np.asarray(p.waveform[0], float)
    x0 = np.asarray(p._period_state[1], float)[:cir.n - 1]
    m = cir.n - 1
    W, c = p._state_event_rows()
    base2, th0 = p._land_fractions(np.diff(ts) / T, [0.706, 0.7107])
    restarts = []

    def FJ(z):
        xx, th = z[:m], z[m:]
        fr, hsens, nodes = p._event_remap(base2, th0, th, T)
        hs_ = fr * T
        tms_ = np.concatenate(([0.0], np.cumsum(hs_)))
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            wk = p._walk('glm', xx, tms_, hs_, T=T, hsens=hsens,
                         capture=set(nodes), keep=True)
        restarts.append(sum(1 for rec in wk.steps if rec.restarted))
        F = np.concatenate((xx - np.asarray(wk.x_end),
                            [float(W[k] @ p._captured[nd][0]) - c[k]
                             for k, nd in enumerate(nodes)]))
        J = np.zeros((m + 2, m + 2))
        J[:m, :m] = np.eye(m) - wk.P
        for k in range(2):
            J[:m, m + k] = -np.asarray(wk.Pk[k]).ravel()
        for k, nd in enumerate(nodes):
            _xj, Pj, Pkj = p._captured[nd]
            J[m + k, :m] = W[k] @ Pj
            for l_ in range(2):
                J[m + k, m + l_] = float(W[k] @ Pkj[l_])
        return F, J
    z = np.concatenate((x0, th0 + np.array([0.001, -0.001])))
    _F0, J0 = FJ(z)
    assert restarts[0] >= 1, restarts
    for i in list(range(m)) + [m, m + 1]:
        zp, zm = z.copy(), z.copy()
        zp[i] += 1e-6
        zm[i] -= 1e-6
        fd = (FJ(zp)[0] - FJ(zm)[0]) / 2e-6
        assert np.linalg.norm(J0[:, i] - fd) < 1e-6 * np.linalg.norm(fd), i

    seed, Tl = _relaxation_oscillator_seed(_comparator_relaxation_oscillator())

    def solve(method, N, se):
        q = PSS(_comparator_relaxation_oscillator(), method=method, reltol=1e-9)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            q.solve(period=Tl, timestep=Tl / N, x0=seed, maxiterations=100,
                    state_events=se)
        assert q.converged, (method, N, se)
        return q
    q3 = solve('glm3', 100, True)
    u3 = solve('glm3', 100, False)
    assert len(q3._state_event_fracs) == 4
    assert abs(q3.period / Tl - 1.0) < 2e-6, q3.period / Tl - 1.0
    assert abs(u3.period / Tl - 1.0) > 1e-3, u3.period / Tl - 1.0
    assert abs(float(np.max(np.abs(q3.floquet_multipliers))) - 1.0) < 1e-3
    q2 = solve('glm2', 200, True)
    assert abs(q2.period / Tl - 1.0) < 2e-4, q2.period / Tl - 1.0


def test_the_landed_edges_drop_the_order_so_a_stiff_state_does_not_ring():
    """A two-step method takes ONE backward-Euler step after each edge the
    grid lands on, as the stepping loop does after a breakpoint.  On a
    STIFF RC (tau = 1e-4 T << h) the resistor current follows the source's
    slope within tau, so a grid point after the ramp starts reads `C/tr`.
    ⚠ Measured before (2026-09-28), when the shooting never dropped: trap
    RANG -- 0.85, 0.73, 0.62 of the current at the first three points after
    an edge, decaying by its stiff factor -- and gear was 0.42 off for a
    step; the uniform grid rings as well (trap 0.83).  With the drop 0.039.
    The drop is keyed to the edge's NODE, which the state-event stage's
    remap moves (`test_the_plain_map_lands_state_events_opened_at_x0`)."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    RR, TAU, TR = 1e3, 1e-4 * T, 0.02 * T
    TD, PW = 0.0125 * T, 0.4 * T
    c = SubCircuit()
    c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=TD, tr=TR, tf=TR, pw=PW,
                     per=T)
    c['R'] = R(1, 2, r=RR)
    c['C'] = C(2, gnd, c=TAU / RR)
    i_ramp = (TAU / RR) / TR
    rows = [str(n) for n in c.nodes]
    for method in ('trap', 'gear'):
        pss = PSS(c, method=method, reltol=1e-10)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / 400, maxiterations=40)
        assert pss.converged and pss.break_events
        ts = np.asarray(pss.waveform[0], float).ravel()
        X = np.asarray(pss.waveform[1], float)
        i = (X[rows.index('1')] - X[rows.index('2')]) / RR
        j0 = int(np.searchsorted(ts, TD * (1 + 1e-9), side='right'))
        err = np.abs(i[j0:j0 + 3] - i_ramp) / i_ramp
        assert np.max(err) < 0.1, (method, err)
    ## `order_drop_at_edges=False` walks through the edges at full order:
    ## trap rings again (the measured 0.85 at the first point)
    pss = PSS(c, method='trap', reltol=1e-10, order_drop_at_edges=False)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, maxiterations=40)
    ts = np.asarray(pss.waveform[0], float).ravel()
    X = np.asarray(pss.waveform[1], float)
    i = (X[rows.index('1')] - X[rows.index('2')]) / RR
    j0 = int(np.searchsorted(ts, TD * (1 + 1e-9), side='right'))
    assert abs(i[j0] - i_ramp) / i_ramp > 0.5, i[j0] / i_ramp
    with pytest.raises(TypeError, match='order_drop_at_edges'):
        PSS(c, method='trap', order_drop_at_edges='yes').solve(
            period=T, timestep=T / 400)


def test_a_staged_solve_that_stalls_falls_back_to_the_one_stage_orbit():
    """`_staged_fallback`: a solve with the state events as Newton unknowns
    that fails solves the one-stage problem from the same seed (its own
    opener) and stages from that orbit.  ⚠ Measured before (2026-09-28):
    trap on the PWM loop from zeros, with the order drop at its landed ramp
    edges, stalled at |F| 3.8 after 1017 evaluations -- the staged solve's
    first stage is the map opened AT `x(0)`, which from that seed does not
    converge, where the one-stage solve with its opener does and the staged
    solve from its orbit converges in 13."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-5
    th_radau = np.array([0.688459, 0.693154, 0.992922, 0.992972])   # radau, 400 pts
    p = PSS(_pwm_loop(T), method='trap', reltol=1e-10)
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        p.solve(period=T, timestep=T / 60, x0=np.zeros(_pwm_loop(T).n - 1),
                maxiterations=100)
    assert p.converged and p.staged_fallback
    assert any('one-stage problem first' in str(r.message) for r in rec)
    assert np.max(np.abs(np.asarray(p._state_event_fracs) - th_radau)) < 1e-3
    ## a solve that stages at once does not take it
    q = PSS(_pwm_loop(T), method='trap', reltol=1e-10,
            order_drop_at_edges=False)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        q.solve(period=T, timestep=T / 60, x0=np.zeros(_pwm_loop(T).n - 1),
                maxiterations=100)
    assert q.converged and not q.staged_fallback
