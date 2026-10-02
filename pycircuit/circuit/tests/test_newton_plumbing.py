"""The Newton's plumbing around a step, trimmed bit for bit (speed round
4, stage C; 2026-10-02).

The reference row went in and out of the state through `np.insert` and
one-element `concatenate`s (9.6 us an iteration in the limiter alone),
the tolerance vectors were built and reduced at every step, the row
names for a failure message too, the Newton's test took every magnitude
twice, the chord took the held Jacobian's every iteration, the history
rings were rebuilt through a one-row array and a view, and both
per-element hooks (`accept_step`, `next_event`) polled every element.
Each item here is the same arithmetic on the same values, or a cache of
a pure function keyed on everything it reads.
"""
import numpy as np
import pytest

from pycircuit.circuit import _tran_history, _tran_newton
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.analysis import insert_row, reduced_row_names
from pycircuit.circuit.circuit import Circuit
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, VSin, gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient

H = 2e-8


def _rect():
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('in', 'a', 'k'):
        c.add_node(n)
    c['vin'] = VSin('in', gnd, v=0.0, va=1.0, freq=1e6)
    c['rs'] = R('in', 'a', r=100.0)
    c['d'] = eh.DiodeSpiceHdl('a', 'k')
    c['rl'] = R('k', gnd, r=1e3)
    c['cl'] = C('k', gnd, c=1e-10)
    return c


def test_a_limited_step_runs_without_the_toolkits_insert(monkeypatch):
    """The limiter's reinserts were `toolkit.insert` (`np.insert`): with
    it made to raise, a limited transient still runs (`uic`: the operating
    point's own limiter keeps `np.insert`, which is DC's business).  Fails
    on the parent, which calls it twice an iteration."""
    def boom(*a, **k):
        raise AssertionError('toolkit.insert called')
    monkeypatch.setattr(numeric, 'insert', boom)
    tr = Transient(_rect(), toolkit=numeric, uic=True, chord_jacobian=False)
    tr.solve(tend=20 * H, timestep=H, fixed_timestep=True)
    assert tr.statistics.newton_iterations > tr.statistics.accepted_steps


@pytest.mark.parametrize('n', [0, 1, 3, 6])
def test_insert_row_is_the_concatenate_form_bit_for_bit(n):
    """Every position, signed zeros and NaN kept, the inserted row +0.0; an
    int vector and a matrix take the toolkit form (float64 result)."""
    rng = np.random.default_rng(n)
    x = rng.standard_normal(6) * 10.0 ** rng.integers(-9, 9, 6)
    x[1] = -0.0
    x[4] = np.nan
    ref = np.concatenate((x[:n], np.array([0.0]), x[n:]))
    got = insert_row(x, n, numeric)
    assert got.tobytes() == ref.tobytes() and got.dtype == np.float64
    assert np.copysign(1.0, got[n]) == 1.0
    xi = np.arange(6)
    assert insert_row(xi, n, numeric).tobytes() == np.concatenate(
        (xi[:n], np.array([0.0]), xi[n:])).tobytes()


def test_the_tolerance_vectors_are_built_once_and_follow_a_parameter(monkeypatch):
    """One `newton_tolerance_vectors` call per run (the parent made two a
    step); a changed `iabstol` is a new key, and the vectors read it."""
    calls = [0]
    real = _tran_newton.newton_tolerance_vectors

    def counting(*a, **k):
        calls[0] += 1
        return real(*a, **k)
    monkeypatch.setattr(_tran_newton, 'newton_tolerance_vectors', counting)
    tr = Transient(_rect(), toolkit=numeric, uic=True)
    tr.solve(tend=20 * H, timestep=H, fixed_timestep=True)
    assert calls[0] == 1
    assert float(tr._newton_abstol_vector()[0]) == float(tr.par.iabstol)
    tr.par.iabstol = 3e-11
    tr.solve(tend=20 * H, timestep=H, fixed_timestep=True)
    assert calls[0] == 2
    assert float(tr._newton_abstol_vector()[0]) == 3e-11
    a, x, ar, xr = tr._newton_tolerances()
    assert ar.shape[0] == a.shape[0] - 1 and xr.shape[0] == x.shape[0] - 1
    assert tr._newton_abstol_vector_reduced() is ar


def test_the_tolerance_vectors_are_read_only():
    """A consumer writing into a cached vector would poison every later
    step: it raises instead (the parent handed out fresh, writable ones)."""
    tr = Transient(_rect(), toolkit=numeric, uic=True)
    tr.solve(tend=5 * H, timestep=H, fixed_timestep=True)
    for v in tr._newton_tolerances():
        with pytest.raises(ValueError):
            v[0] = 1.0


def test_the_reduced_row_names_are_kept_per_circuit_shape():
    tr = Transient(_rect(), toolkit=numeric, uic=True)
    tr.solve(tend=5 * H, timestep=H, fixed_timestep=True)
    names = tr._reduced_row_names()
    assert names is tr._reduced_row_names()
    assert names == reduced_row_names(tr.cir, tr.irefnode)
    assert len(names) == tr.cir.n - 1


def test_the_ring_push_is_the_concatenate_form_bit_for_bit():
    rng = np.random.default_rng(3)
    ring = rng.standard_normal((3, 5))
    v = rng.standard_normal(5)
    v[2] = -0.0
    ref = np.concatenate((np.array([v]), ring))[:-1]
    got = _tran_history._ring_push(v, ring, numeric)
    assert got.tobytes() == ref.tobytes() and got.shape == ring.shape
    assert got.base is None                      # owns its memory
    ## a row of another dtype or shape keeps the toolkit form
    vi = np.arange(5)
    assert _tran_history._ring_push(vi, ring, numeric).tobytes() == \
        np.concatenate((np.array([vi]), ring))[:-1].tobytes()


def _hooked():
    cm.default_toolkit = numeric
    c = SubCircuit()
    for n in ('p', 'a', 'k'):
        c.add_node(n)
    c['vp'] = VPulse('p', gnd, v1=0.0, v2=1.0, td=1e-7, tr=1e-8, tf=1e-8,
                     pw=2e-7, per=5e-7)
    c['r1'] = R('p', 'a', r=1e3)
    c['r2'] = R('a', 'k', r=1e3)
    c['cl'] = C('k', gnd, c=2e-11)
    return c


def test_the_hooks_poll_only_the_elements_that_define_them(monkeypatch):
    """`Circuit`'s no-op `next_event` / `accept_step` are never called
    (the parent polled every element with each); an element that shadows
    a hook on the instance is polled, and its event is the circuit's."""
    seen = {'next_event': 0, 'accept_step': 0}
    base_ne, base_as = Circuit.next_event, Circuit.accept_step

    def ne(self, t):
        seen['next_event'] += 1
        return base_ne(self, t)

    def acc(self, t, x, epar):
        seen['accept_step'] += 1
        return base_as(self, t, x, epar)
    monkeypatch.setattr(Circuit, 'next_event', ne)
    monkeypatch.setattr(Circuit, 'accept_step', acc)
    c = _hooked()
    hits = [0]

    def my_accept(t, x, epar):
        hits[0] += 1
    c['r2'].accept_step = my_accept
    c['r1'].next_event = lambda t: 3.3e-7
    tr = Transient(c, toolkit=numeric, uic=True)
    tr.solve(tend=10 * H, timestep=H, fixed_timestep=True)
    assert seen == {'next_event': 0, 'accept_step': 0}
    assert hits[0] == tr.statistics.accepted_steps + 1    # + the run's start
    ## the full poll, as it was, is the oracle: every element, `max(t, min)`
    for t in (0.0, 2e-7, 3.15e-7, 3.25e-7, 1e-6):
        full = np.maximum(t, min(base_ne(el, t) if type(el).next_event is ne
                                 and 'next_event' not in el.__dict__
                                 else el.next_event(t)
                                 for el in c.elements.values()))
        assert float(c.next_event(t)) == float(full)
    assert float(c.next_event(3.25e-7)) == pytest.approx(3.3e-7)  # r1's shadow
    ## no element declaring an event: `inf`, as the full poll gave
    d = SubCircuit()
    d.add_node('a')
    d['r'] = R('a', gnd, r=1.0)
    assert d.next_event(1.0) == np.inf


def test_a_replaced_element_rebuilds_the_hook_lists():
    """The lists follow the topology: an element replaced by one that
    declares an event is polled after the replacement."""
    c = _hooked()
    assert float(c.next_event(0.0)) == pytest.approx(1e-7)
    c['vp'] = VS('p', gnd, v=1.0)                 # no events any more
    assert c.next_event(0.0) == np.inf
    c['r1'] = VPulse('p', 'a', v1=0.0, v2=1.0, td=4e-7, tr=1e-8, tf=1e-8,
                     pw=2e-7, per=5e-7)
    assert float(c.next_event(0.0)) == float(c['r1'].next_event(0.0)) < np.inf
