"""Tests for DCSweep -- SPICE's .dc (stage 10.1)."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import numeric, gnd, SubCircuit
from pycircuit.circuit.elements import R, VS, Diode
from pycircuit.circuit.dcanalysis import DC, DCSweep


def _divider():
    cir = SubCircuit(toolkit=numeric)
    cir.add_node('a'); cir.add_node('b')
    cir['V1'] = VS('a', gnd, v=0.0)
    cir['R1'] = R('a', 'b', r=1e3)
    cir['R2'] = R('b', gnd, r=1e3)
    return cir


def _diode_chain(n=4):
    cir = SubCircuit(toolkit=numeric)
    names = ['n%d' % i for i in range(n + 1)]
    for nm in names:
        cir.add_node(nm)
    cir['V1'] = VS(names[0], gnd, v=0.0)
    cir['R1'] = R(names[0], names[1], r=1e3)
    for i in range(n - 1):
        cir['D%d' % i] = Diode(names[i + 1], names[i + 2])
    cir['Rl'] = R(names[n], gnd, r=1e6)
    return cir


def test_dcsweep_matches_an_analytic_curve():
    """A resistive divider: the swept answer is exactly half the source."""
    cir = _divider()
    values = np.linspace(0, 5, 11)
    res = DCSweep(cir, toolkit=numeric).solve('V1', 'v', values)
    assert np.allclose(np.asarray(res.sweep_values, dtype=float), values)
    assert np.allclose(np.asarray(res.v('b', gnd), dtype=float), values / 2.0,
                       rtol=1e-9, atol=1e-12)


def test_dcsweep_restores_the_swept_parameter():
    """A sweep must leave the circuit as it found it.

    Otherwise every later analysis silently inherits the last swept value -- the
    same shape of defect stage 8(d) found in TLine, where a DC answer depended on
    whether a transient had run first.
    """
    cir = _divider()
    cir['V1'].ipar.v = 1.25
    DCSweep(cir, toolkit=numeric).solve('V1', 'v', np.linspace(0, 5, 5))
    assert float(cir['V1'].ipar.v) == pytest.approx(1.25)


def test_dcsweep_restores_the_parameter_even_when_a_point_fails():
    """The restore is in a `finally`, and this is what proves it."""
    cir = _divider()
    cir['V1'].ipar.v = 0.75
    with pytest.raises(Exception):
        ## An impossible sweep value type reaches float() and raises mid-sweep.
        DCSweep(cir, toolkit=numeric).solve('V1', 'v', ['not-a-number'])
    assert float(cir['V1'].ipar.v) == pytest.approx(0.75)


def test_dcsweep_rejects_unknown_instance_and_parameter():
    cir = _divider()
    with pytest.raises(ValueError, match='not an instance'):
        DCSweep(cir, toolkit=numeric).solve('NOPE', 'v', [0.0])
    with pytest.raises(ValueError, match='no parameter'):
        DCSweep(cir, toolkit=numeric).solve('V1', 'nope', [0.0])
    with pytest.raises(ValueError, match='empty'):
        DCSweep(cir, toolkit=numeric).solve('V1', 'v', [])


def test_continuation_cuts_the_work_and_not_the_answer():
    """Seeding each point with the previous solution is the point of a sweep.

    Measured on a 4-diode chain over 0..5 V in 101 points: 250 residual
    evaluations with continuation against 1052 without, for the same curve.
    Asserted as a ratio with the answers required to agree, so it says
    "continuation is doing work" rather than pinning today's count.
    """
    import pycircuit.circuit.nrsolver as nrs

    counts = {'n': 0}
    real = nrs.StandardNewton.solve_system

    def patched(self, x0, func, *a, **k):
        seen = {'k': 0}

        def spy(x, *aa, **kk):
            seen['k'] += 1
            return func(x, *aa, **kk)
        out = real(self, x0, spy, *a, **k)
        counts['n'] += seen['k']
        return out

    nrs.StandardNewton.solve_system = patched
    try:
        results, work = {}, {}
        for cont in (True, False):
            counts['n'] = 0
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = DCSweep(_diode_chain(), toolkit=numeric).solve(
                    'V1', 'v', np.linspace(0, 5, 101), continuation=cont)
            results[cont] = np.asarray(res.v('n4', gnd), dtype=float)
            work[cont] = counts['n']
    finally:
        nrs.StandardNewton.solve_system = real

    assert work[True] * 2 < work[False], \
        'continuation saved little: %d evaluations against %d' % (work[True], work[False])
    assert np.abs(results[True] - results[False]).max() < 1e-5, \
        'continuation changed the answer'


def test_dc_x0_defaults_to_the_old_behaviour():
    """`x0=None` must be exactly `zeros`, or DCSweep's seeding would be a change
    in disguise rather than an addition."""
    cir = _diode_chain()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = np.asarray(DC(_diode_chain(), toolkit=numeric).solve().x, dtype=float)
        b = np.asarray(DC(_diode_chain(), toolkit=numeric).solve(
            x0=np.zeros(cir.n)).x, dtype=float)
    assert np.array_equal(a, b)

    with pytest.raises(ValueError, match='unknowns'):
        DC(_diode_chain(), toolkit=numeric).solve(x0=np.zeros(3))


def test_diode_limiter_state_does_not_outlive_its_analysis():
    """Stage 10.1, found while measuring continuation.

    `Diode.limit` stores `_vlim` on the instance and `G` linearises around it, so
    a second DC on the same circuit used to converge in 2 residual evaluations
    where the first took 15 -- a difference that is device state, not the circuit.
    Same defect stage 8(d) found in TLine.history, fixed by the same hook.
    """
    import pycircuit.circuit.nrsolver as nrs
    counts = {'n': 0}
    real = nrs.StandardNewton.solve_system

    def patched(self, x0, func, *a, **k):
        seen = {'k': 0}

        def spy(x, *aa, **kk):
            seen['k'] += 1
            return func(x, *aa, **kk)
        out = real(self, x0, spy, *a, **k)
        counts['n'] += seen['k']
        return out

    nrs.StandardNewton.solve_system = patched
    try:
        cir = _diode_chain()
        runs = []
        for _ in range(3):
            counts['n'] = 0
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                DC(cir, toolkit=numeric).solve()
            runs.append(counts['n'])
    finally:
        nrs.StandardNewton.solve_system = real

    assert len(set(runs)) == 1, \
        'repeated DC on one circuit took %r evaluations -- state is surviving' % runs
