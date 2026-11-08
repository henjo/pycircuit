"""THE RADAU STALL CHECK (`Transient.radau_stall_check`, 2026-11-08;
doc/hdl_roadmap_260824.md, 297; Andreas: "1", from 1 024 unknowns).

An armed Radau attempt whose transform Newton runs `maxiter` iterations
without converging fell back to the coupled Newton, which solves the step
exactly -- and on a 2-D mesh every one was the step over a pulse's edge,
45-63x over its tolerance, rejected by the error test after the coupled
Newton had spent 0.4 s (2 030 unknowns) to 10 s (8 105) on it.  The early
abort cannot see it: at the first iterate the stages have moved as a linear
response, which the embedded pair reads as exact.  At the stalled iterate
the estimate reads 92-130.  So the early abort's check is made again there,
against `radau_stall_check`: over it the attempt is abandoned as an early
abort is, its retry sized by the family's law at the iterate's worst
unknown; else the coupled Newton as before.  'auto', the default, makes it
from `_tran_radau.STALL_CHECK_FROM_N` (1 024) unknowns: on the small gold
decks the coupled Newton costs milliseconds and keeps steps the check turns
away (todd3 1.2-1.7x slower with it on).

What these hold: a stalled attempt over the threshold is abandoned as a
rejection and the coupled Newton not run; a check that never aborts changes
nothing; 'auto' is the check from its size and nothing below it; the C
transform loop makes it as the Python does -- on pattern pages too, where
the loop has no coupled Newton of its own and handed every stall back; the
check rides on the early abort's arm; its threshold is 'auto', 0 or at least
1.  (The mesh here is under the size: the check asked for by number, or the
size lowered.)"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, _tran_radau, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_tcloop_c import _same
from pycircuit.circuit.transient import Transient

#: (300 ns: the pulse's first edge, at 105 ns, and its fall)
_RUN = {'tend': 3e-7, 'timestep': 1e-9}


def _mesh(rows=10, cols=10):
    """Rows x cols inverters, each driven by the node above (row 0 by the
    pulse), every node tied to its four neighbours by 300 k: at the pulse's
    edge the whole first row switches within one step."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    nm = lambda i, j: f'n{i}_{j}'
    for i in range(rows):
        for j in range(cols):
            c.add_node(nm(i, j))
    for i in range(rows):
        for j in range(cols):
            n = nm(i, j)
            c[f'rl{i}_{j}'] = R('vdd', n, r=2e4)
            c[f'M{i}_{j}'] = eh.MosLevel1Hdl(n, 'g0' if i == 0 else nm(i - 1, j), gnd, gnd,
                                             vto=0.5, kp=2e-4)
            c[f'c{i}_{j}'] = C(n, gnd, c=2e-14)
            if j + 1 < cols:
                c[f'rh{i}_{j}'] = R(n, nm(i, j + 1), r=3e5)
            if i + 1 < rows:
                c[f'rv{i}_{j}'] = R(n, nm(i + 1, j), r=3e5)
    return c


def _run(**par):
    """The mesh under Radau's defaults (the C loops where they serve), with
    `par` set: `(the transient, its record, the paths counted)`."""
    circuit.default_toolkit = circuit.numeric
    tr = Transient(_mesh(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
    for k, v in par.items():
        setattr(tr.par, k, v)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**_RUN)
    st = tr.statistics
    rec = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.force_accepts, st.early_aborts)
    return tr, rec, _paths.since(before)


def test_a_stalled_attempt_over_its_tolerance_is_abandoned():
    """The mesh, the check at 2 against the check off: each step over the
    pulse's edge stalls, its estimate there over the threshold -- abandoned
    as a rejection, the coupled Newton not run; with the check off, the
    coupled Newton solved those steps for the error test to reject."""
    tr, _rec, on = _run(radau_stall_check=2)
    off_tr, _r, off = _run(radau_stall_check=0)
    st = tr.statistics
    assert on.get('radau.stall:aborts', 0) == on.get('radau.stall:checks', 0) >= 2, on
    assert not on.get('radau.transform:fallback'), on
    assert off.get('radau.transform:fallback', 0) >= 2 and 'radau.stall:checks' not in off, off
    ## (an abandoned attempt is a rejection, counted where the early abort's
    ## are: the run's statistic holds both)
    assert st.early_aborts == on.get('radau.early:aborts', 0) + on['radau.stall:aborts']
    assert st.rejected_steps >= st.early_aborts
    assert st.force_accepts == off_tr.statistics.force_accepts == 0


def test_a_stall_check_that_never_aborts_changes_nothing():
    """Every stall checked, none over the threshold (infinite: a stalled
    iterate's estimate can be, gm1's is): the run with the check off, byte
    for byte -- the check writes nothing the run reads."""
    _tr, rec, counts = _run(radau_stall_check=float('inf'))
    _o, rec_off, _c = _run(radau_stall_check=0)
    assert counts.get('radau.stall:checks', 0) >= 2 and not counts.get('radau.stall:aborts')
    assert rec == rec_off


def test_auto_is_the_check_from_its_size_and_nothing_below(monkeypatch):
    """'auto' on the mesh's 105 unknowns: no check, the run with it off;
    the size lowered under them: the run with the check at 2."""
    assert Transient(_mesh(), toolkit=circuit.numeric).par.radau_stall_check == 'auto'
    _t, below, c_below = _run()
    _t, off, _c = _run(radau_stall_check=0)
    assert 'radau.stall:checks' not in c_below and below == off
    monkeypatch.setattr(_tran_radau, 'STALL_CHECK_FROM_N', 100)
    _t, above, c_above = _run()
    _t, two, _c = _run(radau_stall_check=2)
    assert c_above.get('radau.stall:aborts', 0) >= 2 and above == two


def test_the_stall_check_rides_on_the_early_abort_s_arm():
    """The early abort off: no attempt armed, no stall checked -- the
    coupled Newton at every stall, as before."""
    _tr, _rec, counts = _run(radau_early_abort=0, radau_stall_check=2)
    assert 'radau.stall:checks' not in counts and counts.get('radau.transform:fallback', 0) >= 2


@pytest.mark.parametrize('pages, lowered', [(False, False), (True, False), (False, True)],
                         ids=['dense', 'pattern-pages', 'auto-from-a-lowered-size'])
def test_the_c_loop_makes_the_stall_check_as_the_python(pages, lowered, monkeypatch):
    """The transform's loop in C and in Python: the same stalls checked and
    abandoned, the same bytes, statistics, counts and state -- made in C
    (no hand-back).  On pattern pages the loop binds no coupled Newton, and
    a stall was a hand-back: the check is what keeps those steps in C.
    'auto' resolved the same in both (`stall_threshold`)."""
    if pages:
        monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    if lowered:
        monkeypatch.setattr(_tran_radau, 'STALL_CHECK_FROM_N', 100)
    ca = _same(_mesh, monkeypatch, steps_min=50, make=None if lowered else {'radau_stall_check': 2},
               **_RUN)
    assert ca.get('radau.stall:aborts', 0) >= 2, ca
    assert not ca.get('tcloop_c:bail:transform') and not ca.get('radau.transform:fallback'), ca


def test_a_stall_checked_and_kept_is_the_python_s(monkeypatch):
    """Every stall checked, none over the threshold (infinite): on to the
    fallback after `tc_stall`'s check -- in C as in the Python, the same
    run."""
    ca = _same(_mesh, monkeypatch, steps_min=50, make={'radau_stall_check': float('inf')}, **_RUN)
    assert ca.get('radau.stall:checks', 0) >= 2 and not ca.get('radau.stall:aborts'), ca
    assert ca.get('radau.transform:fallback', 0) >= 2 and ca.get('tcloop_c:fallbacks', 0) >= 1, ca


@pytest.mark.parametrize('make', [{'radau_stall_check': 0},
                                  {'radau_early_abort': 0, 'radau_stall_check': 2}, None],
                         ids=['threshold-0', 'not-armed', 'auto-below-its-size'])
def test_a_stall_the_c_does_not_check_is_the_python_s(make, monkeypatch):
    """The check off, the attempt not armed, or 'auto' under its size: the
    stall passes `tc_stall` unchecked to the fallback -- in C as in the
    Python, the same run, no stall counted."""
    ca = _same(_mesh, monkeypatch, steps_min=50, make=make, **_RUN)
    assert 'radau.stall:checks' not in ca and ca.get('radau.transform:fallback', 0) >= 2, ca
    ## (the C loop met the stall: its own fallback, or the hand-back)
    assert ca.get('tcloop_c:fallbacks', 0) >= 1, ca


@pytest.mark.parametrize('bad', [0.5, 'on'])
def test_a_threshold_not_auto_0_or_at_least_one_is_refused(bad):
    with pytest.raises(ValueError, match='radau_stall_check'):
        _run(radau_stall_check=bad)
