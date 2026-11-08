"""`save=[...]` (`_tran_result.resolve_save`, `_tran_grid.SavedRows`) and a
transient result's currents (`_tran_result.TransientResult`).

A run that saves keeps the rows named, every accepted point, the bytes of
the run that keeps everything, on the Python's accept and on each C loop's
write-back, with `outputstep` the bytes of its grid; it holds a few points,
not all; a name that does not resolve, or a device whose current the run
cannot give, is refused before the run.  A transient result's currents: a
source's own terminal at each point's time, a diode read at each point (its
current its series resistor's), a device holding charge or a delay line
refused, a plain resistor's as before, bit for bit."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _tran_grid, _tran_result, circuit
from pycircuit.circuit.analysis import CircuitResult
from pycircuit.circuit.elements import C, Diode, ISin, R, SubCircuit, TLine, VSin, gnd
from pycircuit.circuit.tests.test_output_grid import CASES, _solve
from pycircuit.circuit.transient import Transient

#: (the nodes, branches and resistive terminals each case saves)
NAMES = {'_edge': ['d0', 'd3', 'i(vdd.plus)', 'i(rl2.plus)'], '_ramp': ['out', 'Idtmod.idt_node', 'i(vin.plus)', 'i(R1.plus)']}


def _names(case):
    return NAMES[CASES[case][0].__name__]


def _wave(w):
    return np.asarray(getattr(w, 'y', w), float).tobytes()


@pytest.mark.parametrize('case', list(CASES))
def test_the_saved_rows_are_the_whole_run_s(case, monkeypatch):
    """Every saved signal and row the full run's bytes, the times and steps
    the same, through a chunk's edge; the C write-back reached where the
    case runs a C loop."""
    full, tr0, _k = _solve(case, monkeypatch)
    monkeypatch.setattr(_tran_grid, 'SAVE_CHUNK', 7)
    res, tr, k = _solve(case, monkeypatch, save=_names(case))
    assert np.asarray(res.sweep_values).tobytes() == np.asarray(full.sweep_values).tobytes()
    assert tr.statistics.accepted_steps == tr0.statistics.accepted_steps
    assert np.asarray(res.x).tobytes() == np.ascontiguousarray(np.asarray(full.x)[res.rows]).tobytes()
    for name in _names(case):
        get = (lambda r, n=name: r.i(n[2:-1])) if name.startswith('i(') else (lambda r, n=name: r.v(n))
        assert _wave(get(res)) == _wave(get(full)), name
    assert (k > 0) == CASES[case][3]


@pytest.mark.parametrize('case', list(CASES))
def test_saved_on_the_grid_are_the_grid_s(case, monkeypatch):
    """With `outputstep`, the saved rows of the grid stage 1 makes, a block
    every step."""
    step = CASES[case][2]['tend'] / 137
    monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', 1)
    full, _tr, _k = _solve(case, monkeypatch, outputstep=step)
    res, _tr, _k = _solve(case, monkeypatch, outputstep=step, save=_names(case))
    assert np.asarray(res.sweep_values).tobytes() == np.asarray(full.sweep_values).tobytes()
    assert np.asarray(res.x).tobytes() == np.ascontiguousarray(np.asarray(full.x)[res.rows]).tobytes()


def test_the_saved_rows_wait_for_the_shifts(monkeypatch):
    """The idtmod's wraps rewrite the newest three points (its state row is
    saved above, the bytes the full run's); kept one point sooner, the
    store drops a point the rewrite still reaches, and the run fails
    there rather than keeping a value from before it."""
    monkeypatch.setattr(_tran_grid, 'FINISHED_AFTER', _tran_grid.FINISHED_AFTER - 1)
    with pytest.raises(TypeError, match='NoneType'):
        _solve('idtmod, the Python', monkeypatch, save=['Idtmod.idt_node'])


def test_a_saving_run_holds_a_few_points(monkeypatch):
    held = []
    advance = _tran_grid.SavedRows.advance

    def counted(self, step=None):
        advance(self, step)
        held.append(sum(r is not None for r in self.X))
    monkeypatch.setattr(_tran_grid.SavedRows, 'advance', counted)
    _res, tr, _k = _solve('gear2, the C loop', monkeypatch, save=['d0'])
    assert tr.statistics.accepted_steps > 100 and len(held) == tr.statistics.accepted_steps
    assert max(held) <= _tran_grid.FINISHED_AFTER + 1


def _rc_source():
    c = SubCircuit()
    c['I1'] = ISin(gnd, 'a', i=0.0, ia=1e-3, freq=1e3)
    c['R1'] = R('a', gnd, r=1e3)
    c['C1'] = C('a', gnd, c=1e-7)
    return c


def _diode():
    c = SubCircuit()
    c['vs'] = VSin('in', gnd, v=0.0, va=1.0, freq=1e3)
    c['R1'] = R('in', 'a', r=1e3)
    c['D1'] = Diode('a', gnd)
    return c


def _run(c, **par):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(c, toolkit=circuit.numeric, **par)
        return tr.solve(tend=1e-3, timestep=1e-5), tr


@pytest.mark.parametrize('bad, err', [(['nope'], 'not a node'), (['R1.plus'], 'not a node'), (['i(R9.plus)'], 'no terminal'),
                                      (['i(R1.nope)'], 'no terminal'), ('a', 'list of names')])
def test_a_name_save_cannot_keep_is_refused_before_the_run(bad, err, monkeypatch):
    def never(*a, **k):
        raise AssertionError('the run began')
    monkeypatch.setattr(Transient, '_open_run', never)
    with pytest.raises((ValueError, TypeError), match=err):
        _run(_rc_source(), save=bad)


def test_a_source_s_own_terminal_is_its_value_at_each_point():
    """`extract_i` takes the sources at `t=0` unless told: `res.i()` read
    0.0 at every point of a sine current source's terminal."""
    res, _tr = _run(_rc_source())
    t = np.asarray(res.sweep_values, float)
    np.testing.assert_allclose(np.asarray(res.i('I1.plus').y, float), 1e-3 * np.sin(2 * np.pi * 1e3 * t), rtol=1e-12, atol=1e-18)


def test_a_diode_is_read_at_each_point():
    """A stateful limiter's device read at the run's last point read 2e-12 A
    where its series resistor carried 0.42 mA: now KCL holds at every point
    within the Newton's tolerance, and the device's state is put back."""
    c = _diode()
    res, tr = _run(c)
    iR, iD = np.asarray(res.i('R1.plus').y, float), np.asarray(res.i('D1.plus').y, float)
    assert np.abs(iR).max() > 4e-4
    assert np.all(np.abs(iR - iD) <= tr.par.reltol * np.abs(iR) + tr.par.iabstol * 100)
    before = dict(c['D1'].__dict__)
    res.i('D1.plus')
    assert {k: v for k, v in c['D1'].__dict__.items() if 'lim' in k} == {k: v for k, v in before.items() if 'lim' in k}


def test_a_device_holding_charge_is_refused():
    res, _tr = _run(_rc_source())
    with pytest.raises(ValueError, match='holds charge'):
        res.i('C1.plus')


def test_a_delay_line_s_port_is_its_branch_and_what_sums_its_history_is_refused():
    """A `TLine`'s port current is an unknown (a branch row), saved and read
    as one; a current made of its history-dependent stamp -- a subcircuit's
    node it sits on -- is refused (`hidden_state`)."""
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, v=0.0, va=1.0, freq=1e8)
    c['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
    c['RL'] = R('b', gnd, r=50.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        full = Transient(c, toolkit=circuit.numeric).solve(tend=5e-9, timestep=1e-10)
        res = Transient(c, toolkit=circuit.numeric, save=['i(T1.p2)']).solve(tend=5e-9, timestep=1e-10)
    assert _wave(res.i('T1.p2')) == _wave(full.i('T1.p2')) == _wave(CircuitResult.i(full, 'T1.p2'))
    sub = SubCircuit()
    sub['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
    assert 'carries state' in _tran_result._refusal(sub, 'a', 'X1.a')


def test_a_plain_resistor_s_current_is_as_it_was():
    """No sources, no limiter, no charge: `CircuitResult.i`'s bytes."""
    res, _tr = _run(_diode())
    assert _wave(res.i('R1.plus')) == _wave(CircuitResult.i(res, 'R1.plus'))
    assert _wave(res.i('vs.plus')) == _wave(CircuitResult.i(res, 'vs.plus'))


def test_a_signal_not_saved_is_refused_naming_save():
    res, _tr = _run(_diode(), save=['a'])
    assert res.saved == ('a',) and list(res.rows) == sorted(res.rows)
    with pytest.raises(ValueError, match='save'):
        res.v('in')
    with pytest.raises(ValueError, match='save'):
        res.i('R1.plus')


## ---- the currents of devices holding charge, made during the run ----

def _w(r, term):
    return np.asarray(r.i(term).y, float)


@pytest.mark.parametrize('integ, par', [('Gear2Integrator', {}), ('TrapezoidalIntegrator', {}), ('EulerIntegrator', {}),
                                        ('RadauIIA3Integrator', {'radau_transform': False}),
                                        ('RadauIIA3Integrator', {'radau_transform': True})])
def test_a_charge_device_s_current_closes_its_node_s_kcl(integ, par):
    """An RC fed a sine current: the capacitor's current, made from its own
    charge by each step's companion (trap's carries its last current) or,
    under Radau, the step's collocation derivative on its stage charges, and
    the resistor's and the source's add to zero at every point to rounding;
    against the analytic `C dv/dt`, off by the method's own error in `v`
    (over R), no more."""
    from pycircuit.circuit import integrator as I
    res, _tr = _run(_rc_source(), integrator=getattr(I, integ)(), save=['a', 'i(I1.minus)', 'i(R1.plus)', 'i(C1.plus)'], **par)
    t, v = np.asarray(res.sweep_values, float), np.asarray(res.v('a').y, float)
    iI, iR, iC = _w(res, 'I1.minus'), _w(res, 'R1.plus'), _w(res, 'C1.plus')
    assert np.abs(iC).max() > 1e-4
    assert np.abs(iI + iR + iC).max() <= 1e-12 * np.abs(iR).max()
    w, r, c = 2 * np.pi * 1e3, 1e3, 1e-7
    vex = 1e-3 * r / (1 + (w * r * c) ** 2) * (np.sin(w * t) - w * r * c * np.cos(w * t) + w * r * c * np.exp(-t / (r * c)))
    assert np.abs(iC - (1e-3 * np.sin(w * t) - vex / r)).max() <= 1.01 * np.abs(v - vex).max() / r


#: (the MOS chain's node d2 and every device on it, the capacitor and a
#: transistor's other terminals)
CHARGED = ['d2', 'i(rl2.minus)', 'i(M2.d)', 'i(M3.g)', 'i(c2.plus)', 'i(M3.d)', 'i(M0.s)']


@pytest.mark.parametrize('case, unit', [('gear2, the C loop', '_tran_loop_c'), ('radau, the C stage loop', '_tran_stageloop_c'),
                                        ('radau, the C transform loop', '_tran_tcloop_c')])
def test_the_c_loop_and_the_python_make_the_same_currents(case, unit, monkeypatch):
    """The MOS chain with capacitors, each C loop writing back every seven
    steps (the multistep loop its steps' companions, the stage loops their
    stages on the saved devices' rows) against the Python's attempts: every
    current the same bytes; KCL at d2, every device on it saved, within the
    Newton's tolerance."""
    import importlib

    from pycircuit.circuit import _tran_loop_c
    monkeypatch.setattr(_tran_loop_c, 'CAP_MAX', 7)
    on, _tr, k = _solve(case, monkeypatch, save=CHARGED)
    assert k > 0
    monkeypatch.setattr(importlib.import_module('pycircuit.circuit.' + unit), 'ENABLED', False)
    off, tr, k = _solve(case, monkeypatch, save=CHARGED)
    assert k == 0
    for name in CHARGED[1:]:
        assert _w(on, name[2:-1]).tobytes() == _w(off, name[2:-1]).tobytes(), name
    node = sum(_w(on, n[2:-1]) for n in CHARGED[1:5])
    scale = max(np.abs(_w(on, n[2:-1])).max() for n in CHARGED[1:5])
    assert scale > 1e-5 and np.abs(node).max() <= tr.par.reltol * scale + 10 * tr.par.iabstol


def test_saving_currents_changes_nothing_else(monkeypatch):
    full, tr0, _k = _solve('gear2, the C loop', monkeypatch)
    res, tr, _k = _solve('gear2, the C loop', monkeypatch, save=CHARGED)
    assert np.asarray(res.x).tobytes() == np.ascontiguousarray(np.asarray(full.x)[res.rows]).tobytes()
    assert tr.statistics.accepted_steps == tr0.statistics.accepted_steps


@pytest.mark.parametrize('case', ['gear2, the C loop', 'radau, the C stage loop'])
def test_the_currents_on_the_grid_are_their_resample(case, monkeypatch):
    """On a grid a charge-holding device's current is its series at the
    solver's times resampled, bit for bit; a device without charge has its
    current made from the grid's rows, as a full result's (S2a's line)."""
    from pycircuit.circuit.transient import resample_uniform
    step = 4e-8 / 137
    res, _tr, _k = _solve(case, monkeypatch, save=CHARGED)
    monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', 1)
    grid, _tr, _k = _solve(case, monkeypatch, save=CHARGED, outputstep=step)
    assert set(res.currents) == {n[2:-1] for n in CHARGED[2:]}
    for name in CHARGED[2:]:
        _g, want = resample_uniform(res.sweep_values, _w(res, name[2:-1]), step=step)
        assert _w(grid, name[2:-1]).tobytes() == want.tobytes(), name


def test_a_charge_current_this_run_cannot_make_is_refused(monkeypatch):
    """A stage method other than Radau (TR-BDF2's derivative is not a
    collocation's), and a device whose state is rewrapped."""
    from pycircuit.circuit.integrator import TRBDF2Integrator
    with pytest.raises(ValueError, match='not TRBDF2Integrator'):
        _run(_rc_source(), integrator=TRBDF2Integrator(), save=['i(C1.plus)'])
    with pytest.raises(ValueError, match='rewrapped'):
        _solve('idtmod, the Python', monkeypatch, save=['i(Idtmod.iplus)'])


@pytest.mark.parametrize('integ, names, warns', [('RadauIIA3Integrator', ['i(C1.plus)'], True),
                                                 ('RadauIIA3Integrator', ['a', 'i(R1.plus)'], False),
                                                 ('Gear2Integrator', ['i(C1.plus)'], False)])
def test_a_radau_run_saving_charge_currents_says_what_they_are(integ, names, warns):
    """Radau's device currents close a node only as well as its stage solve
    converged: a run saving one says so, once; a multistep run, or a Radau
    run saving no charge current, says nothing of it."""
    from pycircuit.circuit import integrator as I
    from pycircuit.circuit.simwarnings import AccuracyWarning
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        Transient(_rc_source(), toolkit=circuit.numeric, integrator=getattr(I, integ)(), save=names).solve(
            tend=1e-3, timestep=1e-5)
    said = [w for w in W if issubclass(w.category, AccuracyWarning) and 'collocation derivative' in str(w.message)]
    assert len(said) == (1 if warns else 0)
