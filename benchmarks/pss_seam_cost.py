"""What the PSS cold-start seam actually COSTS, as opposed to what it reads.

The LTE report landed 2026-09-01 with three numbers, and the loudest of them
is `max_lte_seam`: on the Q=20 resonator at 100 points/period it reads 38.6x
tolerance for Gear-2 and 15.1x for trapezoidal, against interior per-step
peaks of 0.076 and 0.024.  A ratio of ~500.

⚠ THAT IS AN ESTIMATOR, NOT A COST.  Two facts sit awkwardly beside it:
Gear-2's answer at that grid is only 1.17% low, and it improves 5x when the
grid is halved -- which is not what a defect that dominates the error and
gets WORSE under refinement should do.  So before anyone builds the fix (make
the seeded history part of the shooting unknowns, so the period map is a
fixed point in the state a two-step method needs), the question to answer
is how much of
the answer the seam is actually responsible for.

THE DECOMPOSITION.  `_begin_period` re-seeds a flat charge history at every
shooting iteration -- which is exactly what keeps phi a function of `x0`
alone -- so PSS's answer is the fixed point of a map that cold-starts.  Keep
integrating past the converged solution WITHOUT re-seeding and the ring
buffers fill with the true periodic past; the trajectory then drifts to the
limit cycle the same grid and the same method would produce with no seam.
That gives three numbers on one axis:

    analytic  <---- interior ---->  primed  <-- seam -->  PSS answer

  seam cost      = |primed - PSS answer|
  interior cost  = |primed - analytic|

The circuit is chosen so the third point is exact rather than measured: a
series RLC driven at resonance with a 1 V source has |V_C| = Q volts, so the
analytic peak is 20 V with no reference run to trust.

Anti-vacuity: the primed continuation must actually differ from the cold one,
and the ring buffer must actually hold three distinct charges when it starts
-- both are asserted, because a priming that silently did nothing would show
up as "the seam costs zero", which is the answer this script exists to test.

RESULT (2026-09-01, reltol=1e-9 so the shooting solve's own slack is far
below what is being measured -- at reltol=1e-3 the euler and trap readings
were ~1e-5 V of solve residue, not seam):

    method  npts     PSS     primed   seam cost   interior   seam/interior
    euler    100   8.81540   8.81540    5.1e-12    1.12e+01        0.000
    euler    400  15.20994  15.20994    7.2e-13    4.79e+00        0.000
    gear     100  19.76639  19.89297    1.27e-01   1.07e-01        1.183
    gear     200  19.95451  19.98524    3.07e-02   1.48e-02        2.083
    gear     400  19.99008  19.99735    7.27e-03   2.65e-03        2.742
    trap     100  19.98968  19.98968    1.3e-11    1.03e-02        0.000
    trap     400  19.99957  19.99957    6.3e-13    4.27e-04        0.000

⚠ THE SEAM IS REAL FOR GEAR-2 AND EXACTLY ZERO FOR THE OTHER TWO, and the
report as shipped flagged all three.  `max_lte_seam` read 15.1 for
trapezoidal at npts=100 while its cost is 1.3e-11 V.  The mechanism: only
the third divided difference INSIDE THE ESTIMATOR touches the fabricated
charge for a one-step method, never the companion, so the reading is an
artefact of the measurement and not a property of the answer.  Euler's
companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and `iq_{n-1}`,
and the Euler opening step supplies an `iq` consistent with it (which is
what the first-step order drop is FOR).  Gear-2 reads `q_{n-2}` -- the
entering unknown `x_in` -- and that is the only method that can see it.

WHY GEAR-2 PAYS.  `x_in` is the shooting unknown, and the shooting condition
is `x(0) = x(P)`; it does NOT constrain `x_in` to be the orbit's own
`x(-dt)`.  `x_in` reaches `x(0)` through one order-dropped Euler step, so it
sits O(h^2) off the true history point -- and Gear-2 then reads it AS a
history point.  Confirmed by removing the drop, which makes Gear-2 WORSE at
every grid (err 2.34e-1 -> 4.39e-1 at npts=100, 9.92e-3 -> 2.34e-2 at 400):
the drop is protective, and the 1.27e-01 is the residue it leaves behind.

AND IT DOES NOT GO AWAY.  The seam falls as h^2 (4.13x, 4.23x per halving)
while the interior falls faster (7.25x, 5.57x), so its share of the error
GROWS: 54% of Gear-2's total error at npts=100, 68% at 200, 73% at 400.
That is the case for making the entering history part of the shooting
unknowns -- for multistep methods only.  The prize is measurable: Gear-2's
error at npts=100 would go 2.34e-1 -> 1.07e-1 at no extra cost per
iteration, and the factor grows with refinement.

⚠ ONE EXPERIMENT HERE IS INVALID AND IS KEPT OUT ON PURPOSE.  Re-introducing
a single mechanism (the order drop, or the flat seeding) into the free
continuation does NOT measure that mechanism's cost inside PSS: a shooting
fixed point ABSORBS a once-per-period perturbation by moving `x0`, and a
driven continuation cannot.  Measured, the contradiction is flat -- forcing
the order drop once per period into the trap continuation costs 0.234 V,
while the same drop inside PSS costs 1.3e-11 V.  Only fixed-point against
fixed-point is a fair comparison, which is what `run` does.
"""
import warnings
from copy import copy

import numpy as np

from pycircuit import circuit
from pycircuit.circuit import SubCircuit, gnd, VSin, R, L, C
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.shooting import PSS

## Enough periods for the primed continuation to reach its own limit cycle.
## The deviation decays at the circuit's per-period rate, exp(-pi/Q) = 0.8546,
## so 80 periods take a starting gap to 0.8546**80 = 3.5e-6 of itself.
SETTLE_PERIODS = 80


def q20_rlc(f0=1e3, Q=20.0):
    """Series RLC at resonance: |V_C| = Q * va, exactly."""
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def peak_of(pss, X, node='c'):
    """|V_node| over a collected period, from the reduced state vectors."""
    i = pss.cir.get_node_index(node) - (1 if pss.cir.get_node_index(node)
                                        > pss.irefnode else 0)
    return float(np.max(np.abs(np.asarray([x[i] for x in X], dtype=float))))


def _q_dot(pss, x_full, t):
    """`qdot = -(i(x) + u(t))` -- free from the DAE, exactly, no solve.

    The residual a converged step satisfies IS `i(x) + iq + u = 0`, so the
    companion current is the exact charge derivative at that point.  This is
    the asymmetry the whole seam rests on: the DAE hands over the first
    derivative of the charge and refuses the second (`qddot` needs `xdot`,
    which needs `C^-1`, and C is singular in MNA).  It is also why
    TRAPEZOIDAL has no seam -- it needs only `iq_{-1}`, which is exactly
    this -- while Gear-2 needs a second CHARGE, which no residual equation
    supplies.
    """
    tr = pss._transient()
    i = np.asarray(pss.cir.i(x_full, tr.epar), dtype=float)
    u = np.asarray(pss.cir.u(t, analysis=pss.par.analysis), dtype=float)
    return -(i + u)


def _install_back_extrapolated(pss, x0_full, times, dt):
    """METHOD H: build `q_{-1}` from derivatives instead of solving for it.

        q_{-1} = q_0 - (3h/2) qdot_0 + (h/2) qdot_1     error (5/12) h^3 q3

    Both derivatives are free (`_q_dot`).  `x_1` comes from ONE throwaway
    backward-Euler predictor step; its own O(h^2) error enters with
    coefficient `h/2` and so lands at O(h^3), which is why a first-order
    predictor is enough for a second-order extrapolation.

    Compare the shipped plain path, which is ALGEBRAICALLY the first-order
    member of this family -- its entering charge is exactly
    `q_0 - h qdot_0` (verified to 1.6e-38) -- so this is the same idea
    carried one term further, and the seam should fall as h^3 not h^2.
    """
    tr = pss._transient()
    n = pss.cir.n
    q0 = np.asarray(pss.cir.q(x0_full, tr.epar), dtype=float)
    qd0 = _q_dot(pss, x0_full, times[0])

    ## The predictor: a fresh run seeded at x_0 opens `is_first_step`, so the
    ## step is order-dropped to backward Euler -- exactly what is wanted, and
    ## it needs no history beyond `q_0`.
    tr._begin_run(x0_full, n)
    tr._dt = dt
    x0_red = np.concatenate((np.asarray(x0_full)[:pss.irefnode],
                             np.asarray(x0_full)[pss.irefnode + 1:]))
    x1_red = pss.solve_timestep(x0_red, times[1], dt)
    x1_full = np.asarray(pss._insert_refnode(x1_red), dtype=float)
    qd1 = _q_dot(pss, x1_full, times[1])

    qm1 = q0 - 1.5 * dt * qd0 + 0.5 * dt * qd1

    ## Install: ring index 0 is `q_0` (the previous point for the step to
    ## `h`), index 1 is `q_{-1}`.  Index 2 is never read by a two-step
    ## companion; it is filled to keep the ring well formed.
    tr._begin_run(x0_full, n)
    tr._dt = dt
    tr._qlast = pss.toolkit.array([q0, qm1, qm1][:len(tr._qlast)])
    tr._q_cache = None
    tr._is_first_step = False
    tr._no_history = False
    tr._dt_last = dt
    tr._dt_last2 = None
    return tr


def solve_back_extrapolated(method, npts, f0=1e3, Q=20.0, reltol=1e-9):
    """Shoot on `x_0` alone, with the history back-extrapolated (method H).

    ⚠ THE JACOBIAN HERE IS APPROXIMATE ON PURPOSE.  `d q_{-1}/d x_0` is a
    real term (`C_0 + (3h/2) G_0 - (h/2) G_1 dx_1/dx_0`) and this seeds the
    sensitivity rings the way the plain path does instead -- flat.  That is
    sound for THIS measurement: a Newton's fixed point is set by its
    RESIDUAL, and only its iteration count by its Jacobian.  Shipping the
    formulation would need the exact one; measuring its accuracy does not.
    """
    from pycircuit.circuit import analysis as _an
    period = 1.0 / f0
    circuit.default_toolkit = circuit.numeric
    pss = PSS(q20_rlc(f0, Q), method=method, reltol=reltol)
    ## One ordinary solve first, purely to initialise `irefnode`, the inner
    ## `Transient` and the autonomy verdict.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=period, timestep=period / npts, maxiterations=40)

    times, dt = pss.toolkit.linspace(0.0, period, num=npts, endpoint=True,
                                     retstep=True)
    m = pss.cir.n - 1
    eye = np.asarray(pss.toolkit.eye(m))

    def traverse(x0_red):
        x0_full = np.asarray(pss._insert_refnode(x0_red), dtype=float)
        _install_back_extrapolated(pss, x0_full, times, dt)
        C0 = np.asarray(pss._C_at(x0_red))
        Px, Cs = [eye, eye], [C0, C0]
        Pq = np.zeros((m, m))
        x = copy(x0_red)
        for t in times[1:]:
            x = copy(pss.solve_timestep(x, t, dt))
            Px_new, Pq = pss._step_sensitivity(
                Px, Cs, Pq, np.asarray(pss._Jf), np.asarray(pss._C))
            Px = [Px_new, Px[0]]
            Cs = [copy(np.asarray(pss._C)), Cs[0]]
        return x, Px[0]

    def func(x0_red):
        x_end, Px = traverse(x0_red)
        return np.asarray(x0_red) - np.asarray(x_end), eye - Px

    tol = _an.newton_tolerance_vectors(
        len(pss.cir.nodes), len(pss.cir.branches), pss.par.iabstol,
        pss.par.vabstol, pss.toolkit)[1]
    (tol,) = remove_row_col((tol,), pss.irefnode, pss.toolkit)
    x0_ss, _i, ier, _m = _an.fsolve(
        func, np.zeros(m), maxiter=40, reltol=reltol, abstol=tol, xtol=tol,
        toolkit=pss.toolkit, full_output=True)
    assert ier == 1, 'back-extrapolated %s/%d did not converge' % (method, npts)

    ## The replay must open the way the solve did -- same lesson as the
    ## shipped path, and the reason this is a call rather than a copy.
    x0_full = np.asarray(pss._insert_refnode(x0_ss), dtype=float)
    _install_back_extrapolated(pss, x0_full, times, dt)
    io_ = pss.cir.get_node_index('c')
    io_ -= 1 if io_ > pss.irefnode else 0
    x = copy(x0_ss)
    peak = abs(float(np.asarray(x)[io_]))
    for t in times[1:]:
        x = pss.solve_timestep(x, t, dt)
        peak = max(peak, abs(float(np.asarray(x)[io_])))
    return peak


def solve_shipped(method, npts, f0=1e3, Q=20.0, reltol=1e-9):
    """The shipped formulation, whatever it is for this method.

    For a companion reaching two charges back this solves for `(x_0, x_{-1})`
    jointly; for the others it is the plain path unchanged.  The gate is that
    it REACHES `primed` -- the seam was the only difference, so removing it
    must land on the cycle a real history produces, not merely nearer to it.
    """
    period = 1.0 / f0
    circuit.default_toolkit = circuit.numeric
    pss = PSS(q20_rlc(f0, Q), method=method, reltol=reltol)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=period, timestep=period / npts,
                        maxiterations=40)
    assert pss.converged, 'shipped %s/%d did not converge' % (method, npts)
    return float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel()))), pss


def run(method, npts, f0=1e3, Q=20.0, reltol=1e-3):
    period = 1.0 / f0
    circuit.default_toolkit = circuit.numeric
    pss = PSS(q20_rlc(f0, Q), method=method, reltol=reltol)
    ## THE PLAIN FORMULATION, held here on purpose.  Gear-2 now solves for
    ## its entering history, so measuring the shipped path would measure a
    ## seam of zero and the record of what the seam WAS would be gone.  This
    ## keeps the before, and `solve_shipped` supplies the after.
    pss._solves_history = lambda: False
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=period, timestep=period / npts,
                        maxiterations=40)
    assert pss.converged, '%s/%d did not converge' % (method, npts)
    cold = float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## THE PRIMED CONTINUATION.  `solve` has just finished its replay, so the
    ## `Transient`'s rings hold the real trailing history of the period --
    ## and crucially `_begin_period` is NOT called again, which is the whole
    ## experiment.  Sources are periodic in T, so advancing `t` by whole
    ## periods reproduces the same excitation on the same grid.
    tr = pss._transient()
    qhist = np.asarray(tr._qlast, dtype=float)
    spread = float(np.max(np.abs(qhist[0] - qhist[-1])))
    assert spread > 0.0, \
        'the history is still flat after a replay -- nothing was primed'

    times, dt = pss.toolkit.linspace(0.0, period, num=npts, endpoint=True,
                                     retstep=True)
    x = np.asarray(res['tpss'].x, dtype=float)[:, -1]
    x = np.concatenate((x[:pss.irefnode], x[pss.irefnode + 1:]))
    last = None
    for k in range(SETTLE_PERIODS):
        window = []
        ## ⚠ `times[1:]`, NOT `times` -- npts-1 steps of `P/(npts-1)`, which
        ## is exactly one period.  PSS's own replay takes npts steps because
        ## its FIRST one manufactures `x(0)` from the shooting unknown (which
        ## is the state at `-dt`); the period map proper is the remaining
        ## npts-1.  Repeating that extra step once per period detunes the
        ## drive by 1/(npts-1) -- measured, it took trapezoidal's continued
        ## amplitude to 18.72 V, and `1/sqrt(1 + (2 Q / 99)**2) * 20` is
        ## 18.55, so the artefact was the whole reading.  It is a good
        ## illustration of what this script is guarding against: a
        ## plausible number that measures the harness.
        for t in times[1:]:
            x = pss.solve_timestep(x, t + (k + 1) * period, dt)
            window.append(np.asarray(x, dtype=float))
        last = window
    primed = peak_of(pss, last)

    fixed, apss = solve_shipped(method, npts, f0, Q)
    ## Method H only says anything for a companion that reads two charges;
    ## for the others there is no `q_{-1}` in the recursion to build.
    bx = (solve_back_extrapolated(method, npts, f0, Q)
          if pss._companion_reach() >= 2 else None)
    return dict(method=method, npts=npts, cold=cold, primed=primed,
                analytic=Q, max_lte=pss.max_lte, total_lte=pss.total_lte,
                seam_lte=pss.max_lte_seam, hist_spread=spread,
                fixed=fixed, solved_history=apss.solved_history,
                fixed_seam=apss.max_lte_seam, backextrap=bx)


def main(reltol=1e-9):
    ## reltol=1e-9, not the shipped default: at 1e-3 the shooting solve's own
    ## slack is ~1e-5 V, which is LARGER than the effect being measured for
    ## two of the three methods and reads as a small non-zero seam.
    rows = [run(m, n, reltol=reltol) for m in ('euler', 'gear', 'trap')
            for n in (100, 200, 400)]

    print('Q=20 series RLC at resonance, analytic peak 20 V exactly, '
          'reltol=%g' % reltol)
    print('primed = %d periods continued past the solve with NO re-seed\n'
          % SETTLE_PERIODS)
    hdr = ('%-6s %5s | %9s %9s %9s | %10s %10s | %10s %6s'
           % ('method', 'npts', 'plain', 'primed', 'shipped', 'seam cost',
              'interior', 'err after', 'gain'))
    print(hdr); print('-' * len(hdr))
    for r in rows:
        seam = abs(r['primed'] - r['cold'])
        interior = abs(r['primed'] - r['analytic'])
        before = abs(r['cold'] - r['analytic'])
        after = abs(r['fixed'] - r['analytic'])
        print('%-6s %5d | %9.5f %9.5f %9.5f | %10.3e %10.3e | %10.3e %5.2fx'
              % (r['method'], r['npts'], r['cold'], r['primed'], r['fixed'],
                 seam, interior, after, before / after if after else
                 float('inf')))

    print('\ndoes the fix REACH the cycle a real history produces?')
    for r in rows:
        gap = abs(r['fixed'] - r['primed'])
        print('  %-6s n=%-4d solved_history=%-6s |shipped - primed| = %.2e  %s'
              % (r['method'], r['npts'], r['solved_history'], gap,
                 'seam reported: %s' % r['fixed_seam']))

    ## METHOD H -- the cheap approximate alternative, kept as a measurement.
    ## It shoots on `x_0` ALONE and builds the history from derivatives the
    ## DAE gives away, so it does not double the unknowns.  The question it
    ## answers is whether an approximation can make the seam stop mattering
    ## without the exact formulation's 2m x 2m solve.
    bx = [r for r in rows if r['backextrap'] is not None]
    if bx:
        print('\nmethod H (derivative back-extrapolation, single unknown):')
        hdr = ('  %-6s %5s | %10s %10s | %10s %10s'
               % ('method', 'npts', 'seam plain', 'seam H', 'share plain',
                  'share H'))
        print(hdr); print('  ' + '-' * (len(hdr) - 2))
        for r in bx:
            interior = abs(r['primed'] - r['analytic'])
            sp = abs(r['primed'] - r['cold'])
            sh = abs(r['backextrap'] - r['primed'])
            print('  %-6s %5d | %10.3e %10.3e | %10.3f %10.4f'
                  % (r['method'], r['npts'], sp, sh,
                     sp / interior, sh / interior))
        for a, b in zip(bx, bx[1:]):
            if a['method'] != b['method']:
                continue
            fa = abs(a['backextrap'] - a['primed'])
            fb = abs(b['backextrap'] - b['primed'])
            pa = abs(a['primed'] - a['cold'])
            pb = abs(b['primed'] - b['cold'])
            print('  %s %d -> %d per halving: plain %5.2fx (h^2 = 4), '
                  'H %5.2fx (h^3 = 8)'
                  % (a['method'], a['npts'], b['npts'], pa / pb, fa / fb))

    print('\nwhat the estimator says vs what it costs, where a seam exists'
          ' at all:')
    for r in rows:
        if r['seam_lte'] is None:
            ## Euler and trapezoidal: the companion never reads the entering
            ## point, so the report no longer claims a seam for them.  That
            ## `None` IS the corrected classification, not missing data.
            continue
        seam = abs(r['primed'] - r['cold'])
        interior = abs(r['primed'] - r['analytic'])
        print('  %-6s n=%-4d LTE says seam/interior = %8.1fx, '
              'the answer says %8.4fx'
              % (r['method'], r['npts'], r['seam_lte'] / r['max_lte'],
                 seam / interior if interior else float('inf')))


if __name__ == '__main__':
    main()
