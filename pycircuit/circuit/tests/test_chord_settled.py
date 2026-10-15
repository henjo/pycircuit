"""The chord's settled test (roadmap 136, `nrsolver._chord_settled`, its
C twin newton.c's `chord_settled`): a chord contracting slowly and
MONOTONICALLY is still lambda/(1-lambda) of its increment from the root
when that increment first passes the update test, so it is not accepted
there; an OSCILLATING one is within |lambda|/(1+|lambda|) of it, and is.

A scalar chord with a held slope J0 on f(x) = x - 1 contracts by
lambda = 1 - 1/J0 exactly: J0 = 10 is monotone at 0.9, J0 = 0.625
oscillates at -0.6."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_newton_c, circuit, nrsolver
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import C, R, VPulse
from pycircuit.circuit.nrsolver import ChordNewton, StandardNewton
from pycircuit.circuit.transient import Transient

numeric = circuit.numeric
RELTOL, XTOL = 1e-3, np.array([1e-6])


def _chord(j0):
    """The chord's answer and iterations on f(x) = x - 1 from 0, its slope
    held at j0."""
    def eval_FJ(x):
        return x - 1.0, np.array([[j0]])

    chord = ChordNewton(lambda x: x - 1.0, StandardNewton())
    x, its = chord.solve_system(np.zeros(1), eval_FJ, numeric, RELTOL,
                                np.array([1e-12]), XTOL, 200)
    assert chord.fallbacks == 0
    return float(x[0]), its


def _update_test_alone(j0):
    """Where the update test alone accepts: the iterate and its count."""
    x = 0.0
    for its in range(1, 201):
        dx = -(x - 1.0) / j0
        x_next = x + dx
        tol = RELTOL * max(abs(x_next), abs(x)) + XTOL[0]
        if abs(dx) < tol and abs(x - 1.0) < RELTOL * (abs(j0) * abs(x_next) + abs(x - 1.0)) + 1e-12:
            return x_next, its
        x = x_next
    raise AssertionError('no acceptance')


def test_a_slow_monotone_chord_is_accepted_only_when_settled():
    """lambda 0.9: the update test alone stops ~9 tolerances from the root;
    the settled test runs on to within one."""
    early, early_its = _update_test_alone(10.0)
    assert abs(early - 1.0) > 5 * RELTOL, early
    x, its = _chord(10.0)
    assert abs(x - 1.0) < RELTOL + XTOL[0], x
    assert its > early_its


def test_an_oscillating_chord_is_accepted_where_the_update_test_accepts_it():
    """lambda -0.6: the error after an increment is under half of it, so the
    settled test never refuses -- the same iterate as the update test."""
    early, early_its = _update_test_alone(0.625)
    x, its = _chord(0.625)
    assert its == early_its
    assert x == pytest.approx(early, abs=1e-14)
    assert abs(x - 1.0) < RELTOL


@pytest.mark.parametrize('lam, settled', [
    (0.4, True), (-0.9, True), (0.5, True), (0.55, False), (0.95, False), (1.0, False),
    (1.3, False)])
def test_the_settled_test_reads_the_signed_ratio(lam, settled):
    """At an increment 0.9 of its tolerance: a negative lambda settles, a
    positive one until 0.9 lambda/(1-lambda) reaches 1 (lambda 1/1.9 ~
    0.526), and lambda >= 1 never."""
    tol = np.array([1.0, 1.0])
    prev = np.array([1.0, 0.1])
    now = np.array([lam, 0.0])
    now = now * (0.9 / max(abs(lam), 1e-300))
    prev = prev * (0.9 / abs(lam))
    assert nrsolver._chord_settled(now, tol, prev, tol) is settled


def _falling_clamp():
    """A state-free diode clamp on a falling edge: the chord's held slope
    from the forward-biased seed is stiffer than the root's, so it
    contracts slowly and monotonically (lambda 0.70, 0.81 measured)."""
    c = SubCircuit()
    c['v1'] = VPulse('a', gnd, v1=0.0, v2=5.0, td=1e-7, tr=1e-9, tf=1e-9, pw=2e-7, per=1.0)
    c['r1'] = R('a', 'b', r=1e3)
    c['d1'] = eh.DiodeSpiceHdl('b', gnd, IS=1e-13)
    c['r2'] = R('b', gnd, r=1e3)
    c['cb'] = C('b', gnd, c=1e-12)
    return c


def _run(on, monkeypatch, settled=None):
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', on)
    if settled is not None:
        monkeypatch.setattr(nrsolver, '_chord_settled', settled)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_falling_clamp(), toolkit=numeric, chord_jacobian=True)
        res = tr.solve(tend=4e-7, timestep=2e-9)
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            tr.statistics.newton_iterations, _paths.since(before))


def test_the_c_chord_settles_as_the_python_one(monkeypatch):
    """The C chord (newton.c's `chord_settled`) serves the falling clamp and
    answers byte for byte as the Python; the rule acts there (the Python
    with it off answers otherwise), so the equality sees it."""
    py = _run(False, monkeypatch)
    cc = _run(True, monkeypatch)
    assert cc[3].get('newton_c:served', 0) > 50, cc[3]
    assert (py[0], py[1], py[2]) == (cc[0], cc[1], cc[2])
    off = _run(False, monkeypatch, settled=lambda *a: True)
    assert (off[0], off[1]) != (py[0], py[1])


def test_a_patched_settled_test_takes_the_python_newton(monkeypatch):
    """The C stands in for `_chord_settled` as defined: patched, it declines."""
    got = _run(True, monkeypatch, settled=lambda *a: True)
    assert got[3].get('newton_c:served', 0) == 0
    assert got[3].get('newton_c:patched', 0) > 0
