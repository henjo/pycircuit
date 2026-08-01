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
