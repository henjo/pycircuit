"""Tests for DCSweep -- SPICE's .dc (stage 10.1)."""

import numpy as np
import pytest

from pycircuit.circuit import numeric, gnd, SubCircuit
from pycircuit.circuit.elements import R, VS, Diode
from pycircuit.circuit.dcanalysis import DC, DCSweep
from pycircuit.circuit.tests._warnpolicy import quiet


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


def test_dcsweep_restores_whether_the_swept_parameter_was_given():
    """The sweep's writes GIVE the parameter (a model may branch on
    `$param_given`): one that was not given is left not given, one that was
    is left given."""
    cir = _divider()
    assert not cir['R1'].ipar.is_given('noisy')
    DCSweep(cir, toolkit=numeric).solve('R1', 'r', [1e3, 2e3])
    assert cir['R1'].ipar.is_given('r')
    cir2 = _divider()
    cir2['V1'].ipar.forget('v')
    DCSweep(cir2, toolkit=numeric).solve('V1', 'v', [0.0, 1.0])
    assert not cir2['V1'].ipar.is_given('v')
    assert cir2['V1'].ipar.v == 0.0


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


@pytest.mark.usefixtures('ordinary_dc')
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
            with quiet():
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
    with quiet():
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
            with quiet():
                DC(cir, toolkit=numeric).solve()
            runs.append(counts['n'])
    finally:
        nrs.StandardNewton.solve_system = real

    assert len(set(runs)) == 1, \
        'repeated DC on one circuit took %r evaluations -- state is surviving' % runs


def test_dcsweep_takes_the_dc_tolerances_and_hands_them_to_its_inner_dc():
    """⚠ A gap the `vabstol` default change (039a017) opened and a peer suite
    found: `DCSweep` declared no tolerances and built its inner `DC` with the
    defaults, so `DCSweep(cir, vabstol=1e-12)` raised "parameter vabstol not
    in parameter dictionary" and a sweep could not be asked for a tolerance at
    all -- while the rule for that change was "name the tolerance, do not
    re-pin".  Now: the same Newton parameters as `DC`, same defaults, passed
    through to every point's solve.
    """
    from pycircuit.circuit import dcanalysis as _dcm
    cir = _diode_chain()
    sw = DCSweep(cir, toolkit=numeric, vabstol=1e-12, reltol=1e-6, maxiter=57)
    assert sw.par.vabstol == 1e-12 and sw.par.reltol == 1e-6 and sw.par.maxiter == 57
    ## defaults are DC's defaults, by construction from the same list
    d, s = DC(cir, toolkit=numeric), DCSweep(cir, toolkit=numeric)
    for name in DCSweep.DC_PASSTHROUGH:
        assert getattr(s.par, name) == getattr(d.par, name), name
    ## and the inner DC receives them: spy on the construction
    seen = {}
    real = _dcm.DC

    class Spy(real):
        def __init__(self, *a, **kw):
            seen.update(kw)
            real.__init__(self, *a, **kw)
    _dcm.DC = Spy
    try:
        res = sw.solve('V1', 'v', np.linspace(0.0, 2.0, 3))
    finally:
        _dcm.DC = real
    assert seen['vabstol'] == 1e-12 and seen['reltol'] == 1e-6 \
        and seen['maxiter'] == 57, seen
    assert np.asarray(res.v('n4', gnd), dtype=float).shape[-1] == 3


def _ladder(n_sections=160):
    """A source, a resistor ladder of `n_sections` (over the 150 unknowns
    from which a DC factors on the circuit's structure with KLU) and a diode
    at its end."""
    cir = SubCircuit(toolkit=numeric)
    cir['V1'] = VS('n0', gnd, v=1.0)
    for k in range(n_sections):
        cir[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=10.0)
        cir[f'G{k}'] = R(f'n{k + 1}', gnd, r=1e5)
    cir['D'] = Diode(f'n{n_sections}', gnd)
    return cir


def test_a_dc_per_point_keeps_the_circuit_s_solver():
    """A new DC per point of a sweep -- the usual way to write one -- keeps
    the circuit's default solver: its choice (KLU on the structure, here),
    its analysis and its factorisation; an explicit solver is the DC's
    own; a topology change brings a new one."""
    from pycircuit.circuit.linearsolver import AutoSolver, StructuralKLU
    cir = _ladder()
    a = DC(cir)
    with quiet():
        x = a.solve().x
    solver = a._get_linearsolver()
    assert type(solver) is AutoSolver and isinstance(solver._choice, StructuralKLU)
    cir['V1'].ipar.v = 0.9
    b = DC(cir)
    with quiet():
        b.solve(x0=x)
    assert b._get_linearsolver() is solver
    own = AutoSolver()
    assert DC(cir, linearsolver=own)._get_linearsolver() is own
    cir['Rx'] = R('n1', gnd, r=1e3)
    c = DC(cir)
    with quiet():
        c.solve()
    assert c._get_linearsolver() is not solver


def test_a_kept_klu_solver_answers_as_a_fresh_one_to_the_rounding():
    """A point solved on the kept solver reuses the last point's pivots:
    its answer is a fresh solver's to the rounding."""
    from pycircuit.circuit.linearsolver import AutoSolver
    cir = _ladder()
    with quiet():
        x = DC(cir).solve().x
        cir['V1'].ipar.v = 0.9
        kept = np.asarray(DC(cir).solve(x0=x).x, dtype=float)
        fresh = np.asarray(DC(cir, linearsolver=AutoSolver(
            structure=None)).solve(x0=x).x, dtype=float)
    assert np.allclose(kept, fresh, rtol=1e-12, atol=1e-15)

