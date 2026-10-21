"""SPICE's noncon for the stateful limiter (D-b, 2026-10-15; roadmap 135):
a Newton iteration on which `elements.Diode` LIMITED -- left its junction
at a voltage other than the one its terminals were handed -- is not
converged (`_limiting.NONCON`, read by `StandardNewton`, `DampedNewton`,
`ChordNewton`, `SchurCoupledNewton` and Fang's coupled step).

Before it, a diode clamp with a parallel path came back as the bare
divider: 5 V - 1 k - b, the diode and 1 k from b to ground solved to
2.5 V, every update and residual within tolerance, because `F` read the
diode as the tangent at its clamped `_vlim` (0.24 V), which carries
nothing; a fresh evaluation there reads 9.5e28 A.  The answer is 0.6297 V.

Each answer here is judged by a FRESH circuit's KCL residual at it (no
element state from the solve) or against the SPICE importer's diode, the
same law without state, on the same grid."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _limiting, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.analysis import analysis_kind
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, C, Diode, R, VPulse
from pycircuit.circuit.nrsolver import (
    ChordNewton,
    DampedNewton,
    SchurCoupledNewton,
    StandardNewton,
)
from pycircuit.circuit.transient import Transient

numeric = circuit.numeric


def _clamp(v1=5.0, r2=1e3, kind='stateful', pulsed=False):
    """V1 - 1 k - b, a diode b to ground, R2 b to ground (none: None);
    `pulsed`: V1 steps 0 -> v1 at 1 us, 1 pF on b."""
    c = SubCircuit()
    if pulsed:
        c['v1'] = VPulse('a', gnd, v1=0.0, v2=v1, td=1e-6, tr=1e-9, tf=1e-9,
                         pw=1.0, per=2.0)
    else:
        c['v1'] = VS('a', gnd, v=v1)
    c['r1'] = R('a', 'b', r=1e3)
    c['d1'] = (Diode('b', gnd) if kind == 'stateful'
               else eh.DiodeSpiceHdl('b', gnd, IS=1e-13))
    if r2:
        c['r2'] = R('b', gnd, r=r2)
    if pulsed:
        c['cb'] = C('b', gnd, c=1e-12)
    return c


def _dc(c, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return DC(c, toolkit=numeric, **kw).solve().x


def _fresh_kcl(build, x):
    """max |KCL| of a FRESH circuit from `build()` at `x` (A)."""
    c = build()
    dc = DC(c, toolkit=numeric)
    with analysis_kind(dc.epar, 'dc'):
        f = (np.asarray(c.i(x, dc.epar), dtype=float)
             + np.asarray(c.u(0, dc.epar, analysis='dc'), dtype=float))
    keep = [k for k in range(c.n) if k != dc.irefnode]
    return float(np.max(np.abs(f[keep])))


@pytest.mark.parametrize('solver', [None, DampedNewton()],
                         ids=['standard', 'damped'])
def test_a_clamp_with_a_parallel_path_is_not_the_bare_divider(solver):
    """2.5 V before (the divider, KCL 9.5e28 A there); now the state-free
    diode's answer and a fresh residual at the noise."""
    kw = {} if solver is None else {'nrsolver': solver}
    x = _dc(_clamp(), **kw)
    ref = _dc(_clamp(kind='spice'))
    b = _clamp().get_node_index('b')
    assert abs(x[b] - ref[b]) < 1e-6, (x[b], ref[b])
    assert abs(x[b] - 0.6297) < 1e-4
    assert _fresh_kcl(_clamp, x) < 1e-6


@pytest.mark.parametrize('r2', [None, 1e2, 1e3, 1e4, 1e5, 1e6])
@pytest.mark.parametrize('v1', [0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
def test_every_clamp_of_the_grid_solves_its_circuit(v1, r2):
    """The reach measured before the fix: 14 of these 36 returned the bare
    divider (V1 10 V, R2 100 ohm: 0.909 V where the clamp holds ~0.6)."""
    x = _dc(_clamp(v1, r2))
    assert _fresh_kcl(lambda: _clamp(v1, r2), x) < 1e-6


@pytest.mark.usefixtures('ordinary_dc')
def test_with_the_check_off_the_clamp_is_the_divider_again(monkeypatch):
    """The tests above can fail: without the check the defect is back (DC's
    acceptance check off too, which would refuse it: the next test)."""
    from pycircuit.circuit import pcnr
    monkeypatch.setattr(_limiting, 'NONCON_CHECK', False)
    monkeypatch.setattr(pcnr, 'DC_ACCEPT_RESIDUAL', float('inf'))
    x = _dc(_clamp())
    assert abs(x[_clamp().get_node_index('b')] - 2.5) < 1e-3
    assert _fresh_kcl(_clamp, x) > 1e20


@pytest.mark.usefixtures('ordinary_dc')
def test_with_the_check_off_dcs_acceptance_check_refuses_the_clamp(monkeypatch):
    """The clamp passes the Newton test but not DC's acceptance check
    (`pcnr.dc_residual`: the devices at the node voltages, a 1e28 A residual
    at b): the solve raises where it returned the divider."""
    from pycircuit.circuit.analysis import NoConvergenceError
    monkeypatch.setattr(_limiting, 'NONCON_CHECK', False)
    with pytest.raises(NoConvergenceError, match='does not solve the circuit'):
        _dc(_clamp())


def _tran(c, check=True, monkeypatch=None, tend=3e-6, **kw):
    """v(b) on the fixed 2 ns grid to `tend`: (t, v)."""
    if not check:
        monkeypatch.setattr(_limiting, 'NONCON_CHECK', False)
    skw = kw.pop('solve', {})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(c, toolkit=numeric, **kw).solve(
            tend=tend, timestep=2e-9, fixed_timestep=True, **skw)
    return (np.asarray(res.sweep_values, dtype=float),
            np.asarray(res.v('b'), dtype=float))


def _off(got, ref):
    """max |v - v_ref| at the reference's grid times (a run that cut a
    step holds them all, its extra points between them)."""
    (t, v), (tr, vr) = got, ref
    return float(np.max(np.abs(np.interp(tr, t, v) - vr)))


CASES = [(5.0, 1e3), (10.0, 1e2), (20.0, 1e4)]
PATHS = {'newton': {}, 'fang': {'solve': {'coupled_lte': True}}}


@pytest.mark.parametrize('path', sorted(PATHS))
@pytest.mark.parametrize('v1, r2', CASES)
def test_a_fixed_step_transient_through_the_clamp_is_the_state_free_diodes(v1, r2, path):
    """The clamp driven 0 -> V1 on one grid: before the fix the stateful
    diode passed false points at the edge, 0.25 .. 11.8 V off; now within
    2e-8 V of the state-free diode."""
    got = _tran(_clamp(v1, r2, pulsed=True), **dict(PATHS[path]))
    ref = _tran(_clamp(v1, r2, kind='spice', pulsed=True))
    assert _off(got, ref) < 1e-6


@pytest.mark.parametrize('v1, r2', CASES)
def test_the_chord_through_the_clamp_passes_no_false_point(v1, r2):
    """`chord_jacobian`'s Newton (`ChordNewton`): before the fix 0.25 ..
    11.9 V off (one run ENDED at 9.84 V); after it (roadmap 135) still
    1.8 .. 2.4 mV, the chord accepting an iterate whose increments
    contracted monotonically at 0.97 -- ~30 tolerances from the root.  With
    the chord's settled test (roadmap 136, `test_chord_settled`) it is the
    state-free diode's chord to 1e-7 V (both 5e-6 V or less from the full
    Newton: the chord's own tolerance)."""
    got = _tran(_clamp(v1, r2, pulsed=True), chord_jacobian=True)
    ref = _tran(_clamp(v1, r2, kind='spice', pulsed=True), chord_jacobian=True)
    assert _off(got, ref) < 1e-6


@pytest.mark.parametrize('kw', [{}, {'chord_jacobian': True},
                                {'solve': {'coupled_lte': True}}],
                         ids=['newton', 'chord', 'fang'])
def test_with_the_check_off_the_transient_passes_false_points(kw, monkeypatch):
    """The transient tests above can fail: without the check every path
    passes a false point at the 5 V case's edge (1.37 V off).  (To 1.2 us:
    the 20 V case, 11.8 V off, cuts steps for 83 s without the check.)"""
    got = _tran(_clamp(5.0, 1e3, pulsed=True), check=False,
                monkeypatch=monkeypatch, tend=1.2e-6, **dict(kw))
    ref = _tran(_clamp(5.0, 1e3, kind='spice', pulsed=True), tend=1.2e-6)
    assert _off(got, ref) > 1.0


def _bumping_limiter(k):
    """A limiter that moves nothing and reports a limited stateful device on
    its first `k` calls."""
    calls = [0]

    def limiter(x_new, x_old):
        calls[0] += 1
        if calls[0] <= k:
            _limiting.NONCON[0] += 1
        return x_new
    return limiter


def _solve_linear(solver, limiter):
    """x - 1 = 0, from 0: converged after one update unless refused."""
    def eval_FJ(x):
        return x - 1.0, np.eye(1)
    if solver == 'schur':
        def eval_FJh(x, h):
            return (np.array([x[0] - 1.0]), np.eye(1), np.zeros(1),
                    h - 1.0, np.zeros(1), 1.0)
        (x, _h), its = SchurCoupledNewton().solve_system(
            (np.zeros(1), 1.0), eval_FJh, numeric, 1e-6, 1e-12, 1e-12, 20,
            limiter=limiter)
        return x, its
    s = {'standard': StandardNewton(), 'damped': DampedNewton(),
         'chord': ChordNewton(lambda x: x - 1.0, StandardNewton())}[solver]
    return s.solve_system(np.zeros(1), eval_FJ, numeric, 1e-6, 1e-12, 1e-12,
                          20, limiter=limiter)


@pytest.mark.parametrize('solver', ['standard', 'damped', 'chord', 'schur'])
def test_every_newton_refuses_an_iteration_its_limiter_limited(solver, monkeypatch):
    """Three limited iterations: no answer before the fourth -- and with the
    check off the first converged iterate is taken."""
    x, its = _solve_linear(solver, _bumping_limiter(3))
    assert abs(x[0] - 1.0) < 1e-12 and its >= 4, its
    monkeypatch.setattr(_limiting, 'NONCON_CHECK', False)
    _x, its_off = _solve_linear(solver, _bumping_limiter(3))
    assert its_off < 4, its_off


def test_a_state_free_limiter_never_counts():
    """The SPICE importer's diode moves the node itself (its update is the
    limited one): it never touches the count -- nor do the gold decks,
    which have no `elements.Diode`."""
    n0 = _limiting.NONCON[0]
    _dc(_clamp(kind='spice'))
    assert _limiting.NONCON[0] == n0


class _Autodiff:
    """`numeric`, saying it differentiates (as `JAXToolkit` does)."""

    def __getattr__(self, name):
        return getattr(numeric, name)

    def supports(self, capability):
        return capability == 'autodiff'


def test_a_diode_counts_only_where_it_reads_its_clamp():
    """A limited junction counts on a numeric toolkit; under an autodiff
    toolkit `i` / `G` read the node voltage, `_vlim` is never read, and
    nothing counts.  An unlimited step (2 VT or less) counts nowhere."""
    epar = DC(_clamp(), toolkit=numeric).epar
    for toolkit, counted in ((numeric, 1), (_Autodiff(), 0)):
        d = Diode('b', gnd)
        d.toolkit = toolkit
        n0 = _limiting.NONCON[0]
        d.limit(np.array([0.3, 0.0]), np.array([0.3, 0.0]), epar)
        assert _limiting.NONCON[0] == n0
        d.limit(np.array([2.5, 0.0]), np.array([0.3, 0.0]), epar)
        assert _limiting.NONCON[0] - n0 == counted
        assert d._vlim < 1.0
