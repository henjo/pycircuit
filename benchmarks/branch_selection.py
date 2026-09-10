"""One netlist, one grid, one tolerance -- THREE different answers.

WHAT THIS IS.  At a critical point of a DAE the solution can be genuinely
non-unique (Lamour/März/Tischendorf Thm 3.53: "there are TWO solutions passing
through").  The `im D` falsifier next door measured a HARMLESS critical point,
where the order collapses but uniqueness survives.  This one is not harmless,
and it takes TWO ingredients rather than one:

  * `rank C` DROPS -- here `C = c0 V^2`, which vanishes at `V = 0`;
  * the equilibrium there is REPELLING -- here a NEGATIVE conductance.

With a passive conductance the equilibrium attracts, the field is one-sided
Lipschitz and forward uniqueness is safe however badly `C` degenerates.  A
negative conductance is not an exotic ingredient: it is what the active device
in an oscillator supplies.

⚠⚠ AND THE KNOB IS NOT THE GRID.  The analytic non-uniqueness becomes
MULTIPLICITY OF ROOTS OF THE STEP EQUATION -- implicit Euler from `v_prev = 0`
gives `z(c0 z^2 / 3h + g) = 0`, which has THREE roots when `g < 0`, namely `0`
and `+-sqrt(-3hg/c0)` -- and the nonlinear solver picks one SILENTLY.  Every
grid picks the same one, so a grid-refinement or grid-offset probe returns a
clean, convincing, WRONG null.  The knob that exposes it is the SOLVER'S
INITIAL GUESS.  (Construction relayed from docs-46, measured here.)

MEASURED, `V` at `t = 1` starting from `V = 0`, by the seed given to the first
step's Newton:

    N       seed -1     seed -0.1    seed 0    seed +0.1   seed +1
    200    -1.420240   -1.420240    0.000000   1.420240   1.420240
    800    -1.416027   -1.416027    0.000000   1.416027   1.416027
    3200   -1.414744   -1.414744    0.000000   1.414744   1.414744

The spread is 2.83 and does NOT shrink -- 2.8405 / 2.8321 / 2.8295.  The
non-trivial branches converge to `+-sqrt(2)`, which is exact: `w' = (3w)^(1/3)`
integrates to `w = ((2/3) 3^(1/3) t)^(3/2)`, so `V(1) = (3w)^(1/3) = sqrt(2) =
1.414214` against a measured 1.414744 at N = 3200.  With `g = +1` the same
sweep returns identically 0.000000 at every seed and every refinement.

⚠ A REPELLING EQUILIBRIUM CANNOT BE ARRIVED AT.  An earlier version of this
file drove an orbit "through" the degenerate point with a current source and
found the predictor made no difference -- on orbits that NEVER CROSSED ZERO.
`min|V|` equalled `|v0|` in every run: with `g < 0` the origin repels, so a
forward trajectory can only LEAVE it.  The branch point is reachable only by
STARTING there, which is the configuration above.  That was this measurement's
instrument failure, and the tell was asking whether the fixture engaged the
condition at all rather than trusting a null.

DOES THE STAGE PREDICTOR (shipped 2026-09-10) CHANGE THE CHOICE?  No -- gear,
TR-BDF2, radau and GLM3 all take the trivial branch with it on and with it
off, because the first step from `V = 0` has no history and the predictor
declines.  That is worth knowing rather than assuming: the predictor changed
every Newton seed in the tree, and the seed is exactly what selects the branch.
"""
import warnings

import numpy as np

import pycircuit.circuit.circuit as _cc
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import G
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter, ddt)
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (Gear2Integrator, TRBDF2Integrator,
                                          RadauIIA3Integrator, GLM3Integrator)

warnings.simplefilter('ignore')
_cc.default_toolkit = numeric


class CubicCap(Behavioural):
    """`q = c0 V^3/3`, so `C = c0 V^2` vanishes at `V = 0`."""
    instparams = [Parameter(name='c0', desc='coefficient', unit='F',
                            default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, ddt(c0 * b.V ** 3 / 3)),)      # noqa: F821


def build(g):
    """`g < 0` makes `V = 0` REPELLING -- the second ingredient."""
    c = SubCircuit()
    c.add_node('a')
    c['cq'] = CubicCap('a', gnd, c0=1.0)
    c['g'] = G('a', gnd, g=g)
    return c


def march(g, npts, seed=None, predictor='on', cls=Gear2Integrator, tend=1.0):
    """March from `V = 0`.  `seed`, if given, is the NEWTON'S initial guess on
    the FIRST step -- the knob that selects the branch.  The grid is untouched.
    """
    cir = build(g)
    ia = cir.get_node_index('a')
    h = tend / npts
    tr = Transient(cir, integrator=cls(), reltol=1e-12)
    tr.irefnode = cir.get_node_index(gnd)
    x = np.zeros(cir.n)
    tr.epar.t = 0.0
    tr._begin_run(x, cir.n)
    state = {'j': 0}
    orig = Transient._predict_state

    def seeded(self, ttarget, extra=(), deg=None):
        if state['j'] == 1 and seed is not None:
            v = np.zeros(cir.n)
            v[ia] = seed
            return v
        return orig(self, ttarget, extra, deg)
    Transient._predict_state = seeded
    Transient.stage_predictor = predictor
    try:
        for j in range(1, npts + 1):
            state['j'] = j
            tr._dt_last = tr._dt if j > 1 else None
            tr._dt = h
            tr.epar.t = j * h
            x, _f, _J, _ = tr.solve_timestep(x, j * h)
            tr._push_history(x)
    finally:
        Transient._predict_state = orig
        Transient.stage_predictor = 'on'
    return float(np.asarray(x, dtype=float)[ia])


SEEDS = (-1.0, -0.1, 0.0, 0.1, 1.0)


def branch_sweep():
    print('=== V(t=1) from V=0, by the FIRST STEP\'S NEWTON SEED ===')
    print('    exact non-trivial branch: V(1) = sqrt(2) = 1.414214')
    for g, lbl in ((-1.0, 'g = -1  (V=0 REPELLING)'),
                   (1.0, 'g = +1  (V=0 ATTRACTING -- the control)')):
        print('  %s' % lbl)
        print('  %6s %s' % ('N', '   '.join('seed=%+.2f' % s for s in SEEDS)))
        for npts in (200, 800, 3200):
            v = [march(g, npts, s) for s in SEEDS]
            print('  %6d %s   spread %.4f'
                  % (npts, '   '.join('%9.6f' % x for x in v),
                     max(v) - min(v)))
        print()


def predictor_changes_nothing():
    """⚠ The predictor changed every Newton seed in this tree, and the seed is
    what selects the branch.  Measured rather than assumed."""
    print('=== does the stage predictor change which branch is taken? ===')
    print('%-8s %6s %16s %16s %12s'
          % ('method', 'N', 'predictor OFF', 'predictor ON', 'same?'))
    for cls, name in ((Gear2Integrator, 'gear'), (TRBDF2Integrator, 'trbdf2'),
                      (RadauIIA3Integrator, 'radau'),
                      (GLM3Integrator, 'glm3')):
        for npts in (200, 800):
            a = march(-1.0, npts, None, 'off', cls)
            b = march(-1.0, npts, None, 'on', cls)
            print('%-8s %6d %16.9f %16.9f %12s'
                  % (name, npts, a, b, 'yes' if abs(a - b) < 1e-9 else 'NO'))
    print()


def the_orbit_must_engage():
    """⚠⚠ THE INSTRUMENT CHECK THAT CAUGHT THIS FILE'S OWN FAILURE.  Driving an
    orbit "through" the degenerate point does not work: with `g < 0` the origin
    REPELS, so `min|V|` never falls below `|v0|` and the crossing never
    happens.  A null measured there says nothing."""
    from pycircuit.circuit.elements import IS
    print('=== can an orbit be driven THROUGH the branch point? ===')
    print('%8s %8s %12s %12s %10s' % ('v0', 'i_src', 'min|V|', 'V end',
                                      'crossed?'))
    for v0, iamp in ((0.6, -0.35), (0.6, -1.5), (0.3, -1.5), (0.1, -2.0)):
        cir = build(-1.0)
        cir['is'] = IS('a', gnd, i=iamp)
        ia = cir.get_node_index('a')
        tr = Transient(cir, integrator=Gear2Integrator(), reltol=1e-12)
        tr.irefnode = cir.get_node_index(gnd)
        x = np.zeros(cir.n)
        x[ia] = v0
        tr.epar.t = 0.0
        tr._begin_run(x, cir.n)
        traj = [v0]
        for j in range(1, 401):
            tr._dt_last = tr._dt if j > 1 else None
            tr._dt = 1.0 / 400
            tr.epar.t = j / 400.0
            x, _f, _J, _ = tr.solve_timestep(x, j / 400.0)
            tr._push_history(x)
            traj.append(float(np.asarray(x, dtype=float)[ia]))
        t = np.array(traj)
        print('%8.2f %8.2f %12.3e %12.6f %10s'
              % (v0, iamp, float(np.min(np.abs(t))), t[-1],
                 'YES' if np.any(np.sign(t[:-1]) != np.sign(t[1:])) else 'no'))
    print()


if __name__ == '__main__':
    the_orbit_must_engage()
    branch_sweep()
    predictor_changes_nothing()
