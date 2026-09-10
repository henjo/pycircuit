# -*- coding: utf-8 -*-
"""The stage predictor: one construction, every integrator family.

Every implicit integrator here used to start its Newton from a value at the
WRONG TIME.  Measured on the fixture below, as a fraction of one step's own
state motion: coupled Radau IIA(3) 1.00 (`x_n` for all three stages at once),
ESDIRK43 0.50 and TR-BDF2 0.59 (the previous stage), Gear-2 and trapezoidal
1.00 (`x_n`, one solve per step).  `Transient._predict_state` replaces all of
them with the polynomial through the recorded nodes nearest the target time.

⚠⚠ THE INSTRUMENT IS THE HARD PART.  `Diode` CANNOT measure a stage predictor:
its `G` linearises around a STORED `_vlim` (its own docstring, stage 13-2), so
its Newton is seed-blind -- measured, an ALL-ZEROS seed gives the iteration
histogram [436, 244, 108], to the count, that the exact seed gives.  A linear
circuit is worse still, forcing a one-step Newton.  `_expg_fixture` is a
state-free exponential with no limiting, which is what it takes.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, R, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import transient as TR
from pycircuit.circuit import nrsolver as NR
from pycircuit.circuit.integrator import (RadauIIA3Integrator,
                                          ESDIRK43Integrator,
                                          TRBDF2Integrator, Gear2Integrator,
                                          TrapezoidalIntegrator,
                                          GLM3Integrator, GLM4Integrator)

PER = 1e-3


def _expg_fixture(per=PER, va=0.8):
    """A driven RC with a STATE-FREE exponential conductance -- `i` and `G`
    both functions of the passed `x`, no stored state, no limiting, so the
    Newton's iteration count responds to its seed."""
    import sympy
    import pycircuit.circuit.circuit as _cc
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter)

    class ExpG(Behavioural):
        instparams = [Parameter(name='IS', desc='sat', unit='A',
                                default=1e-12),
                      Parameter(name='VT', desc='thermal', unit='V',
                                default=0.026)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, IS * (sympy.exp(b.V / VT) - 1)),)  # noqa: F821

    _cc.default_toolkit = numeric
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=va, freq=1.0 / per)
    c['rs'] = R('a', 'b', r=50.0)
    c['nl'] = ExpG('b', gnd, IS=1e-12, VT=0.026)
    c['cl'] = C('b', gnd, c=1e-9)
    c['rl'] = R('b', gnd, r=1e4)
    return c


def _counted(cls, mode, npts, va=0.8, fixed=True, reltol=1e-9):
    """One transient; returns (device `i` evaluations, worst solve's Newton
    iterations, waveform time, waveform value)."""
    log = []
    solve_orig = NR.StandardNewton.solve_system

    def wrapped(self, x0, eval_FJ, *a, **k):
        x, it = solve_orig(self, x0, eval_FJ, *a, **k)
        log.append(it)
        return x, it
    NR.StandardNewton.solve_system = wrapped
    TR.Transient.stage_predictor = mode
    n = {'i': 0}
    cir = _expg_fixture(PER, va)
    fi = cir.i
    cir.i = lambda *a, **k: (n.__setitem__('i', n['i'] + 1), fi(*a, **k))[1]
    try:
        tr = Transient(cir, integrator=cls(), reltol=reltol)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tr.solve(refnode=gnd, tend=PER, timestep=PER / npts,
                           fixed_timestep=fixed)
    finally:
        NR.StandardNewton.solve_system = solve_orig
        TR.Transient.stage_predictor = 'on'
    wv = res.v('b')
    return (n['i'], int(max(log)) if log else 0,
            np.asarray(wv.x[0], dtype=float), np.asarray(wv.y, dtype=float))


## the measured reduction per family, at 200 points and 0.8 V; the gate is set
## well inside it so a modest regression still fails
FAMILIES = [
    (RadauIIA3Integrator, 'radau', 0.08),
    (ESDIRK43Integrator, 'esdirk43', 0.08),
    (TRBDF2Integrator, 'trbdf2', 0.05),
    (Gear2Integrator, 'gear2', 0.10),
    (TrapezoidalIntegrator, 'trap', 0.10),
    (GLM3Integrator, 'glm3', 0.08),
    (GLM4Integrator, 'glm4', 0.15),
]


@pytest.mark.parametrize('cls,name,gain', FAMILIES)
def test_the_stage_predictor_cuts_every_familys_newton_work(cls, name, gain):
    """EVERY family, not just the one it was written for: the sequential-stage
    methods (ESDIRK43, TR-BDF2, the GLMs), the COUPLED one (Radau IIA(3), whose
    three stages are one Newton and so has no stage of its own to read), and
    the multistep ones (Gear-2, trapezoidal, which have no stages at all and
    extrapolate their accepted history instead).

    The `'off'` control is the pre-predictor seed and runs in the same test, so
    this measures the change and not the machine.  On a FIXED grid both runs
    share their time points exactly, so the answer is compared bit for bit
    rather than to a tolerance -- a seed change must not select a different
    root, and that is a real risk on a multi-rooted stage system.
    """
    i_off, it_off, t_off, y_off = _counted(cls, 'off', 200)
    i_on, it_on, t_on, y_on = _counted(cls, 'on', 200)
    assert i_on < (1.0 - gain) * i_off, (name, i_on, i_off)
    assert t_on.shape == t_off.shape and np.allclose(t_on, t_off)
    assert np.max(np.abs(y_on - y_off)) < 1e-11, \
        (name, float(np.max(np.abs(y_on - y_off))))
    assert it_on <= it_off, (name, it_on, it_off)


@pytest.mark.parametrize('cls,name', [(RadauIIA3Integrator, 'radau'),
                                      (ESDIRK43Integrator, 'esdirk43'),
                                      (Gear2Integrator, 'gear2')])
def test_the_stage_predictor_also_pays_on_the_adaptive_path(cls, name):
    """⚠ THE DEFAULT PATH IS THE ADAPTIVE ONE, and it reaches the predictor by
    a different door: the stage-method adaptive driver deliberately calls no
    `_push_history` (a stage method reads no charge rings), so it accepts a
    step without touching any of the other history this class keeps.  A
    predictor hung off `_push_history` alone fires ZERO times there.  ⚠ It is
    RADAU that catches that, not every method here: a sequential-stage method
    still has its own converged stages to predict from within the step, so
    ESDIRK43 keeps most of its gain with the history dead.  The COUPLED and
    MULTISTEP paths have no such fallback and go fully inert.

    Also asserts what an adaptive run may NOT do: the step count must not grow.
    A better seed changes the converged value at the last digits, that moves
    the LTE estimate, and a controller that then took MORE steps would have
    spent the saving.
    """
    i_off, _, t_off, _ = _counted(cls, 'off', 200, fixed=False)
    i_on, _, t_on, _ = _counted(cls, 'on', 200, fixed=False)
    assert i_on < 0.92 * i_off, (name, i_on, i_off)
    assert len(t_on) <= len(t_off) + 2, (name, len(t_on), len(t_off))


@pytest.mark.parametrize('cls,name', [(Gear2Integrator, 'gear2'),
                                      (ESDIRK43Integrator, 'esdirk43'),
                                      (RadauIIA3Integrator, 'radau'),
                                      (GLM4Integrator, 'glm4')])
def test_the_stage_predictor_actually_fires(cls, name):
    """⚠⚠ THE GATE THAT CATCHES A SILENTLY INERT PREDICTOR, and it caught two
    real defects that every cost measurement above would have reported only as
    "no gain":

    * the MULTISTEP path recorded no node of its own -- the stage methods get
      theirs next to `_rk_Y`, and without its own line the history never
      reached two entries: 200 calls, 0 predictions;
    * DUPLICATE NODE TIMES made the fit singular, and they are the normal case
      -- a stiffly accurate method's last stage IS the step it ends, and an
      ESDIRK's explicit first stage IS the state it starts from -- so 37% of
      ESDIRK43's stages silently kept the old seed.

    Both raise nothing and change no answer.  Only a firing rate sees them.
    """
    orig = TR.Transient._predict_state
    st = {'call': 0, 'hit': 0}

    def wrapped(self, ttarget, extra=(), deg=None):
        st['call'] += 1
        r = orig(self, ttarget, extra, deg)
        if r is not None:
            st['hit'] += 1
        return r
    TR.Transient._predict_state = wrapped
    TR.Transient.stage_predictor = 'on'
    try:
        tr = Transient(_expg_fixture(), integrator=cls(), reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr.solve(refnode=gnd, tend=PER, timestep=PER / 200,
                     fixed_timestep=True)
    finally:
        TR.Transient._predict_state = orig
    assert st['call'] > 100, (name, st['call'])
    assert st['hit'] > 0.9 * st['call'], (name, st['hit'], st['call'])


def test_the_stage_predictor_is_exact_on_a_polynomial():
    """STRUCTURAL, no circuit and no integrator: the predictor INTERPOLATES the
    nodes it is given, so on a trajectory that IS a polynomial of its degree it
    returns the value exactly.  That is what says the weights are a Lagrange
    evaluation at the target time and not a fit with a constant in it.

    ⚠ Asserted on a NON-UNIFORM node set, because keying the fit on absolute
    times is the whole reason one predictor serves an adaptive multistep method
    and a fixed-step GLM alike -- a uniform set would pass on a formulation
    that silently assumed constant spacing.
    """
    tr = Transient(_expg_fixture(), integrator=RadauIIA3Integrator())
    times = [0.0, 1.1e-6, 1.7e-6, 3.0e-6, 3.4e-6, 5.0e-6]
    for deg in (1, 2, 3, 4):
        for k in range(deg + 1):
            def traj(tt, _k=k):
                return np.array([tt ** _k, 1.0 - 2.0 * tt ** _k])
            tr._pred_reset()
            for tt in times:
                tr._pred_hist.append((tt, traj(tt)))
            ## strictly inside the node span, where no clamp can bind, so this
            ## measures the weights alone
            for target in (2.0e-6, 2.6e-6, 4.1e-6):
                got = tr._predict_state(target, deg=deg)
                want = traj(target)
                assert np.max(np.abs(got - want)) < 1e-9 * max(
                    float(np.max(np.abs(want))), 1e-12), (deg, k, target)


def test_the_stage_predictor_clamps_to_the_linear_extrapolation():
    """⚠⚠ THE PROPERTY THAT DECIDES WHETHER THIS IS A SPEED-UP OR A
    REGRESSION.  A polynomial continued past its last node can leave the region
    the circuit visits, and on an exponential device a 3x overshoot is
    `exp(3 dV / VT)`: measured, Gear-2 on a coarse grid cost +17% device
    evaluations and a worst solve of 27 Newton iterations against the old
    seed's 8, unclamped.

    The bound is the LINEAR prediction -- from the newest node, moving at the
    rate the last step moved -- times `PRED_CLAMP`.  Gated on a trajectory that
    ACCELERATES hard, where a cubic shoots far past, and on one that does not,
    where the clamp must not bind at all.
    """
    tr = Transient(_expg_fixture(), integrator=RadauIIA3Integrator())
    ts = [0.0, 1e-6, 2e-6, 3e-6]

    tr._pred_reset()
    for tt, v in zip(ts, [0.0, 1.0, 3.0, 9.0]):        # accelerating
        tr._pred_hist.append((tt, np.array([v])))
    motion = 9.0 - 3.0
    ahead = 3.0                                        # 3 steps past the last
    got = float(tr._predict_state(6e-6, deg=3)[0])
    assert got <= 9.0 + tr.PRED_CLAMP * motion * ahead + 1e-9, got
    ## the unclamped cubic is far above that, or this gate proves nothing
    tv = np.array(ts)
    raw = float(np.polyval(np.polyfit(tv, [0.0, 1.0, 3.0, 9.0], 3), 6e-6))
    assert raw > 9.0 + tr.PRED_CLAMP * motion * ahead, raw

    tr._pred_reset()
    for tt, v in zip(ts, [0.0, 1.0, 2.0, 3.0]):        # straight line
        tr._pred_hist.append((tt, np.array([v])))
    assert abs(float(tr._predict_state(4e-6, deg=3)[0]) - 4.0) < 1e-9


def test_a_rejected_steps_stages_never_become_predictor_nodes():
    """A rejected step's stages are samples of a trajectory the run then threw
    away, and they sit BEYOND the step that is eventually accepted (the driver
    halves `dt` and retries), so a predictor that read them would extrapolate
    from points ahead of where it is.  They are promoted at the ACCEPT sites
    only, and this asserts the consequence: no node is ever later than the last
    accepted time.

    Run with a tolerance tight enough to force rejections, and the run is
    checked to HAVE rejected steps -- otherwise this passes on a run that never
    exercised the case.
    """
    tr = Transient(_expg_fixture(va=2.0), integrator=TRBDF2Integrator(),
                   reltol=1e-10)
    seen = []
    orig = Transient._pred_note

    def noted(self, t, x, stages=()):
        orig(self, t, x, stages)
        seen.append((t, [float(e[0]) for e in self._pred_hist]))
    Transient._pred_note = noted
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr.solve(refnode=gnd, tend=PER / 4, timestep=PER / 200,
                     fixed_timestep=False)
    finally:
        Transient._pred_note = orig
    assert tr.statistics.rejected_steps > 0, 'no rejection was exercised'
    for taccept, nodes in seen:
        assert max(nodes) <= taccept + 1e-15, (taccept, max(nodes))


def test_the_stage_predictor_leaves_a_wrapping_state_alone():
    """⚠⚠ A PERIODIC ROW IS NOT A TRAJECTORY A POLYNOMIAL CAN FIT.  It folds by
    its modulus, and that fold is a DISCONTINUITY in exactly the curve being
    fitted -- the gauge shift keeps the recorded nodes in one gauge, but the
    fold can also fall between the newest node and the target, and then the fit
    runs straight across it.

    ⚠ MEASURED as a convergence failure, not a wobble: `Idtmod` with the wrap
    landing exactly ON a grid point, where the period map is genuinely
    discontinuous and the PSS is documented to converge anyway
    (`test_a_state_reset_needs_no_saltation_but_grid_alignment_is_a_cliff`),
    stopped converging at all once those rows were predicted.  They keep the
    old seed; every other row still gets the prediction, which this also
    asserts, or "leave it alone" would be indistinguishable from "give up".
    """
    from pycircuit.circuit.elements import VS, Idtmod
    from pycircuit.circuit.shooting import PSS
    from pycircuit.circuit import circuit as _c

    _c.default_toolkit = _c.numeric
    per = 1e-3

    def build(ic):
        c = SubCircuit()
        for nn in ('in', 'out', 'f'):
            c.add_node(nn)
        c['vin'] = VS('in', gnd, v=2000.0)
        c['I'] = Idtmod('in', gnd, 'out', gnd, modulus=1.0, ic=ic)
        c['Rf'] = R('out', 'f', r=1e3)
        c['Cf'] = C('f', gnd, c=1e-7)
        c['Rl'] = R('out', gnd, r=1e5)
        return c

    ## the grid-ALIGNED wrap, which is the case that failed
    TR.Transient.stage_predictor = 'on'
    p = PSS(build(0.0), method='trap', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p.solve(period=per, timestep=per / 1200, maxiterations=60)
    assert p.converged

    ## and the rule itself: the periodic row keeps the newest node, the others
    ## are predicted
    tr = Transient(_expg_fixture(), integrator=RadauIIA3Integrator())
    tr._pred_reset()
    for k, tt in enumerate([0.0, 1e-6, 2e-6]):
        tr._pred_hist.append((tt, np.array([float(k), 10.0 * k])))
    tr._periodic_rows = []
    free = tr._predict_state(3e-6, deg=2)
    assert abs(free[0] - 3.0) < 1e-9 and abs(free[1] - 30.0) < 1e-9, free
    tr._periodic_rows = [(1, 1.0, 0.0)]
    held = tr._predict_state(3e-6, deg=2)
    assert abs(held[0] - 3.0) < 1e-9, held      # untouched
    assert held[1] == 20.0, held                # the newest node, not 30.0
