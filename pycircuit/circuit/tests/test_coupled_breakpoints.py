"""STAGE 12B, GATE 12-4 -- the coupled path must honour breakpoints.

`_solve_coupled` had no breakpoint handling at all: no `next_event`, no
truncation, no order drop.  That is worse on this path than it would be on the
standard one, because Fang's Figure 3 has no rejection branch -- a coupled step
that runs past a pulse edge cannot back up from it.

It also silently broke the invariant `p` rests on.  `TimeFunction.dfdt` returns
the RIGHT-hand limit at a corner, which is exact only because a step always
starts on the corner rather than straddling it -- and that is true only if the
solver truncates to breakpoints.  On the coupled path it did not, so `p` could
be handed a derivative from the wrong segment.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VPulse
from pycircuit.circuit.transient import Transient

TD, TR, PW, TF, PER = 1e-5, 1e-7, 2e-5, 1e-7, 5e-5
TEND = 1.2e-4


def _pulsed_rc():
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=TD, tr=TR, tf=TF, pw=PW,
                     per=PER)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def _edges():
    out = []
    for k in range(3):
        b = k * PER
        out += [b + TD, b + TD + TR, b + TD + TR + PW, b + TD + TR + PW + TF]
    return [e for e in out if e < TEND]


def _run(coupled):
    tran = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=TEND, timestep=1e-6, coupled_lte=coupled)
    return tran, np.asarray(res.v('b').x, dtype=float).ravel()


@pytest.mark.parametrize('coupled', [False, True], ids=['standard', 'coupled'])
def test_every_pulse_edge_lands_on_a_time_point(coupled):
    """Both paths must place a time point exactly on each edge."""
    _tran, t = _run(coupled)
    for e in _edges():
        miss = float(np.min(np.abs(t - e)))
        assert miss < 1e-12, \
            'edge at %g s missed by %g s (%s path)' \
            % (e, miss, 'coupled' if coupled else 'standard')


def test_the_coupled_path_does_not_solve_away_its_own_truncation():
    """The regression that made the first breakpoint fix pointless.

    `fang_timestep` solves for `h`, so on a truncated step it happily replaced
    the step size that had been chosen to land on the edge and walked off it
    again: 0 of 10 edges hit, worst miss 1.24e-7 s -- the entire rise time.  A
    truncated step's size is IMPOSED, so the LTE equation is dropped for it and
    only the circuit is solved.

    Asserted through the public path rather than by inspecting `hold_h`, so it
    still holds if the mechanism is renamed.
    """
    _tran, t = _run(True)
    worst = max(float(np.min(np.abs(t - e))) for e in _edges())
    assert worst < 1e-12, 'worst edge miss %g s' % worst


def test_both_paths_agree_on_the_waveform_through_the_edges():
    """A step that straddles an edge integrates through a discontinuity.

    Compared on the MEDIAN, not the maximum. The two paths choose different
    grids, so comparing them requires interpolating one onto the other, and
    across a 100 ns rise that interpolation error dominates everything else --
    measured at 4.1e-2 on the maximum against 9.8e-4 on the median. The maximum
    here would be measuring `np.interp`, not the solver.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r_std = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5).solve(
            tend=TEND, timestep=1e-6, coupled_lte=False)
        r_cpl = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5).solve(
            tend=TEND, timestep=1e-6, coupled_lte=True)

    ts = np.asarray(r_std.v('b').x, dtype=float).ravel()
    vs = np.asarray(r_std.v('b').y, dtype=float).ravel()
    tc = np.asarray(r_cpl.v('b').x, dtype=float).ravel()
    vc = np.asarray(r_cpl.v('b').y, dtype=float).ravel()

    d = np.abs(np.interp(ts, tc, vc) - vs)
    assert float(np.median(d)) < 5e-3, \
        'paths disagree, median %g V' % float(np.median(d))


def test_the_coupled_path_records_the_breakpoints_it_hits():
    """`breakpoints_hit` was silently zero on this path.

    A statistic that is always zero is worse than one that is absent: it reads
    as evidence that nothing happened. Order drops are NOT asserted here --
    the default integrator is Euler, which has no lower order to drop to, so
    that counter is legitimately zero on both paths for this circuit.
    """
    tran, _t = _run(True)
    assert tran.statistics.breakpoints_hit >= len(_edges()), \
        'recorded %d breakpoints across %d pulse edges' \
        % (tran.statistics.breakpoints_hit, len(_edges()))


## ---------------------------------------------------------------------------
## GATE 12-4, second input: `fixed_timestep` on the coupled path.

def test_fixed_timestep_keeps_the_grid_on_the_coupled_path():
    """The caller's grid wins over the method's step-size solve.

    `fixed_timestep` and Fang's method are not really in conflict -- one exists
    to CHOOSE the step size and the other says the caller already has. So the
    grid is kept and the LTE equation is dropped on every step, exactly as for a
    breakpoint-truncated one: the circuit is still solved coupled, it just has
    nothing left to solve for.

    Silently adapting anyway would return output points the caller did not ask
    for, which is the defect stage 4h fixed on the standard path.
    """
    step = 1e-6
    tend = 2e-5
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5).solve(
            tend=tend, timestep=step, coupled_lte=True, fixed_timestep=True)

    t = np.asarray(res.v('b').x, dtype=float).ravel()
    dt = np.diff(t)
    assert np.allclose(dt, step, rtol=1e-9), \
        'grid was not uniform: steps ranged %g .. %g' % (dt.min(), dt.max())


def test_fixed_timestep_and_adaptive_differ_on_the_coupled_path():
    """Guard against the test above passing because nothing adapts anyway."""
    step = 1e-6
    tend = 2e-5
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        adaptive = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5).solve(
            tend=tend, timestep=step, coupled_lte=True)
    dt = np.diff(np.asarray(adaptive.v('b').x, dtype=float).ravel())
    assert not np.allclose(dt, step, rtol=1e-9), \
        'the adaptive coupled path produced a uniform grid, so the fixed-step ' \
        'test above proves nothing'


## ---------------------------------------------------------------------------
## GATE 12-4, fourth input: a caller-injected step controller.

def _rc():
    from pycircuit.circuit.elements import VSin
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def test_an_injected_solution_controller_is_the_one_used():
    """Not silently replaced by the path's own.

    Before this, `tran.step_controller = X` looked honoured on the coupled path
    and did nothing -- the same class of defect as a documented feature that does
    not exist.
    """
    from pycircuit.circuit.stepcontroller import SolutionLTEController
    tran = Transient(_rc(), toolkit=numeric, reltol=1e-5)
    ctrl = SolutionLTEController()
    tran.step_controller = ctrl
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=2e-4, timestep=1e-5, coupled_lte=True)
    ## The private fallback must never have been built.
    assert getattr(tran, '_fang_controller', None) is None


def test_an_incompatible_controller_is_refused_with_the_reason():
    """On this path the step law is Fang's, so another controller's cannot run.

    Accepting one and quietly using it only for `relref` would be the same silent
    no-op in a different costume, so it raises and says why.
    """
    from pycircuit.circuit.stepcontroller import IntegralController
    tran = Transient(_rc(), toolkit=numeric, reltol=1e-5)
    tran.step_controller = IntegralController()
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=2e-4, timestep=1e-5, coupled_lte=True)
    assert 'IntegralController' in str(exc.value)
    assert 'coupled_lte=False' in str(exc.value)


def test_the_standard_paths_own_controller_is_not_mistaken_for_an_injection():
    """THE REGRESSION THAT BROKE ELEVEN TESTS.

    `_solve` auto-creates `self.step_controller = IntegralController()` when the
    caller supplied none. Any object that ran the adaptive path first then
    presents that to the coupled path -- and the first version of the check above
    refused it, rejecting a controller nobody had asked for.
    """
    tran = Transient(_rc(), toolkit=numeric, reltol=1e-5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=2e-4, timestep=1e-5, coupled_lte=False)
    assert tran.step_controller is not None      # auto-created by _solve
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=2e-4, timestep=1e-5, coupled_lte=True)
    assert len(np.asarray(res.v('b').x, dtype=float).ravel()) > 10
