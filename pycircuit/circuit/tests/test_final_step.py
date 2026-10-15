"""No rounding-level final step (roadmap 137).  (1) The run's end merges a
source corner within `minbreak`: `_StepEvents._next_breakpoint(t, tend)`
returns `tend` for a corner that close, as SPICE keeps its stop time in the
breakpoint table; loop.c's `next_breakpoint` is its twin.  (2) An adaptive
step ending within `minbreak` short of the end ends at it, and a residue
under `minbreak` left at the end ends the run (`_SteppingLoop.where`, loop.c's
`where`) -- the adaptive twin of the fixed grid's residue rule.

Before, an end one ulp past a corner (a `.tran` stop of 100n read as
100 x 1e-9 against a pulse's 100E-9: hussamp, nand) made the loop land on the
corner and then step 1.3e-23 s; capped steps summed two ulps short of the end
(schmitecl) left a 2.1e-22 s step.  Such a step's companion conductances are
~1e23 and its matrix singular: the Newton returned noise, converged by luck,
and with the standing gmin failed for 100 iterations before a rescue."""
from types import SimpleNamespace

import numpy as np

from pycircuit.circuit import circuit
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.tests.test_loop_c import _pulsed, _run, _same
from pycircuit.circuit.transient import Transient


def _corner(after):
    """The pulse's corner after `after`, as the circuit names it."""
    return float(_pulsed().next_event(after))


def test_the_corner_and_the_end_are_one_ulp_apart_here():
    """The fixture is what it claims: a corner, and an end one ulp past it."""
    c = _corner(1.05e-8)
    tend = float(np.nextafter(c, np.inf))
    assert 0.0 < tend - c < 1e-22 and 1.0e-8 < c < 1.2e-8


def test_an_end_one_ulp_past_a_corner_is_landed_on_in_one_step(monkeypatch):
    """The step that would land on the corner lands on `tend`: the run ends
    there, and no step is a rounding residue."""
    tend = float(np.nextafter(_corner(1.05e-8), np.inf))
    got, _c = _run(_pulsed, False, monkeypatch, tend=tend, timestep=1e-9, fixed_timestep=False)
    t = np.frombuffer(got['t'], dtype=float)
    assert t[-1] == tend
    assert np.min(np.diff(t)) > 1e-18, np.min(np.diff(t))


def test_an_end_one_ulp_before_a_corner_ends_the_run_at_the_end(monkeypatch):
    """The other side: the corner just past `tend` -- the run ends at `tend`."""
    tend = float(np.nextafter(_corner(1.05e-8), -np.inf))
    got, _c = _run(_pulsed, False, monkeypatch, tend=tend, timestep=1e-9, fixed_timestep=False)
    t = np.frombuffer(got['t'], dtype=float)
    assert t[-1] == tend
    assert np.min(np.diff(t)) > 1e-18


def test_the_c_loop_ends_the_run_as_the_python_one(monkeypatch):
    """loop.c's `next_breakpoint` merges the end as `_next_breakpoint` does:
    the same run, byte for byte, the C loop serving."""
    tend = float(np.nextafter(_corner(1.05e-8), np.inf))
    _same(_pulsed, monkeypatch, steps_min=5, tend=tend, timestep=1e-9, fixed_timestep=False)


def test_next_breakpoint_merges_only_a_corner_within_minbreak_of_the_end():
    """`tend` returned for a corner within `minbreak` of it, either side; a
    corner farther off, or no `tend`, as before."""
    tr = Transient(_pulsed(), toolkit=circuit.numeric)
    c = _corner(1.05e-8)
    t = 1.05e-8
    assert tr._next_breakpoint(t) == c
    assert tr._next_breakpoint(t, float(np.nextafter(c, np.inf))) == float(np.nextafter(c, np.inf))
    assert tr._next_breakpoint(t, float(np.nextafter(c, -np.inf))) == float(np.nextafter(c, -np.inf))
    far = c + 1e3 * tr.par.minbreak
    assert tr._next_breakpoint(t, far) == c


def _loop(t, h, tend, fixed=False):
    """Just what `_SteppingLoop.where` reads: no breakpoint ahead."""
    tr = SimpleNamespace(par=SimpleNamespace(minbreak=1e-14), _is_first_step=False)
    family = SimpleNamespace(next_breakpoint=lambda _t: float('inf'))
    return SimpleNamespace(tr=tr, t=t, h=h, tend=tend, fixed=fixed, timestep=h, landing=False,
                           order_drop=False, imposed=False, family=family)


def test_a_step_ending_an_ulp_or_two_short_of_the_end_ends_at_it():
    """schmitecl: a corner landed on, then 5 ns steps (the cap) summed to
    4.999999999999999e-07 against an end of 500 x 1e-9 -- a 2.1e-22 s step
    left; now the step that falls short by under `minbreak` ends at `tend`."""
    t, h = 4.95e-07, 5e-9
    tend = float(np.nextafter(np.nextafter(t + h, np.inf), np.inf))
    lp = _loop(t, h, tend)
    assert trmod._SteppingLoop.where(lp) is True
    assert lp.h == tend - t and t + lp.h == tend and lp.clamped


def test_a_residue_under_minbreak_at_the_end_ends_the_run():
    """Whatever leaves the run an ulp short of its end, the run is over."""
    tend = 5.000000000000001e-07
    lp = _loop(float(np.nextafter(tend, -np.inf)), 5e-9, tend)
    assert trmod._SteppingLoop.where(lp) is False


def test_a_step_farther_from_the_end_and_a_fixed_grid_are_as_before():
    """Only within `minbreak`: a step ending 1e-13 s short keeps its size,
    and a fixed grid keeps its own residue rule."""
    t, h = 4.95e-07, 5e-9
    lp = _loop(t, h, t + h + 1e-13)
    assert trmod._SteppingLoop.where(lp) is True and lp.h == h
    tend = float(np.nextafter(np.nextafter(t + h, np.inf), np.inf))
    lp = _loop(t, h, tend, fixed=True)
    assert trmod._SteppingLoop.where(lp) is True and lp.h == h
