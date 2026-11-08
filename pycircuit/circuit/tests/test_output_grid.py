"""`outputstep`'s grid made as the run goes (`_tran_grid.OutputGrid`): the
result is `resample_uniform` of the whole run's points, every byte, on the
Python's accept and on each C loop's write-back, with a block the run's every
step and with one block at its end; the periodic shifts' rewrites seen (a lag
one point short is another result); the run holding a few points, not all;
and on hand-made runs the edges: a run ending short of `tend`, past it, a
rounding below a grid point, too short to interpolate."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _tran_grid, _tran_loop_c, circuit
from pycircuit.circuit.elements import VS, Idtmod, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    Gear2Integrator,
    RadauIIA3Integrator,
)
from pycircuit.circuit.tests.test_loop_c import _mos
from pycircuit.circuit.transient import Transient, resample_uniform


def _edge():
    c = _mos(6, cap=1e-13)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=3e-9, tr=1e-10, tf=1e-10, pw=6e-9, per=2e-8)
    return c


def _ramp():
    """1 V into an Idtmod of modulus 1: its state wraps every second, and
    each wrap rewrites the run's last three points."""
    c = SubCircuit()
    c['vin'] = VS('in', gnd, v=1.0)
    c['R1'] = R('out', gnd, r=1e3)
    c['Idtmod'] = Idtmod('in', gnd, 'out', gnd, modulus=1.0, offset=-0.5)
    return c


#: (build, Transient's settings, solve's, the C write-back it must reach)
CASES = {
    'gear2, the C loop': (_edge, lambda: {'integrator': Gear2Integrator()}, {'tend': 4e-8, 'timestep': 2.5e-10}, True),
    'radau, the C stage loop': (_edge, lambda: {'integrator': RadauIIA3Integrator(), 'radau_transform': False},
                                {'tend': 4e-8, 'timestep': 2.5e-10}, True),
    'radau, the C transform loop': (_edge, lambda: {'integrator': RadauIIA3Integrator(), 'radau_transform': True},
                                    {'tend': 4e-8, 'timestep': 2.5e-10}, True),
    'idtmod, the Python': (_ramp, lambda: {'integrator': EulerIntegrator(), 'uic': True},
                           {'tend': 5.0, 'timestep': 1e-2, 'fixed_timestep': True}, False),
}


def _solve(case, monkeypatch, **par):
    """`(result, statistics, the C write-backs' points)`."""
    build, make, run, c_loop = CASES[case]
    if c_loop and _tran_loop_c.driver() is None:
        pytest.skip('the stepping loop unit did not build')
    seen = []
    append = _tran_loop_c._append_points

    def counted(loop, Xo, To, k, n):
        seen.append(k)
        return append(loop, Xo, To, k, n)
    with monkeypatch.context() as mp:
        mp.setattr(_tran_loop_c, '_append_points', counted)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr = Transient(build(), toolkit=circuit.numeric, **make(), **par)
            res = tr.solve(**run)
    return res, tr, sum(seen)


def _whole_run(case, monkeypatch, step):
    """`resample_uniform` of the run's own points, as `outputstep` was made
    when the run held them all."""
    res, tr, _k = _solve(case, monkeypatch)
    return resample_uniform(res.sweep_values, res.x, step=step), tr


def _bytes(t, x):
    return np.asarray(t, float).tobytes(), np.asarray(x, float).tobytes()


@pytest.mark.parametrize('per', [137, 2003])
@pytest.mark.parametrize('case', list(CASES))
def test_the_grid_made_as_the_run_goes_is_the_whole_run_s_resample(case, per, monkeypatch):
    """A block every step (`GRID_BLOCK_BYTES` 1) and one block at the end,
    on a grid coarser than the steps and one finer: the bytes of
    `resample_uniform` on every point the run made."""
    step = CASES[case][2]['tend'] / per
    (grid, want), tr0 = _whole_run(case, monkeypatch, step)
    for block in (1, _tran_grid.GRID_BLOCK_BYTES):
        monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', block)
        res, tr, k = _solve(case, monkeypatch, outputstep=step)
        assert _bytes(res.sweep_values, res.x) == _bytes(grid, want), (case, block)
        assert tr.statistics.accepted_steps == tr0.statistics.accepted_steps
        assert (k > 0) == CASES[case][3], f'{case}: {k} points through a C write-back'


def test_the_run_holds_a_few_points_not_all(monkeypatch):
    """With a block a step, the points the run holds stay a handful: the
    finished ones handed to the grid and dropped."""
    held = []
    advance = _tran_grid.OutputGrid.advance

    def counted(self):
        advance(self)
        held.append(sum(r is not None for r in self.X))
    monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', 1)
    monkeypatch.setattr(_tran_grid.OutputGrid, 'advance', counted)
    _res, tr, _k = _solve('gear2, the C loop', monkeypatch, outputstep=4e-8 / 137)
    assert tr.statistics.accepted_steps > 100 and len(held) == tr.statistics.accepted_steps
    assert max(held) <= _tran_grid.FINISHED_AFTER + 3


def test_the_shifts_rewrites_are_waited_for(monkeypatch):
    """The idtmod's wraps rewrite the newest three points after their step
    is accepted: a grid made one point sooner reads them before the rewrite
    and is another result -- so the case above sees the lag."""
    step = 5.0 / 2003
    (grid, want), tr0 = _whole_run('idtmod, the Python', monkeypatch, step)
    assert tr0._periodic_rows
    monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', 1)
    monkeypatch.setattr(_tran_grid, 'FINISHED_AFTER', _tran_grid.FINISHED_AFTER - 1)
    res, _tr, _k = _solve('idtmod, the Python', monkeypatch, outputstep=step)
    assert _bytes(res.sweep_values, res.x)[0] == _bytes(grid, want)[0]
    assert _bytes(res.sweep_values, res.x)[1] != _bytes(grid, want)[1]


def _hand_made(T, step, tend, block, monkeypatch, n=3, seed=0):
    """A run of the given times (after 0) fed point by point: the grid's
    `(times, values)` and the whole run's `resample_uniform`."""
    rng = np.random.default_rng(seed)
    rows = [rng.standard_normal(n) for _ in range(len(T) + 1)]
    X, TT = [rows[0].copy()], []
    monkeypatch.setattr(_tran_grid, 'GRID_BLOCK_BYTES', block)
    g = _tran_grid.OutputGrid(X, TT, step, tend)
    for r, t in zip(rows[1:], T):
        X.append(r.copy())
        TT.append(float(t))
        g.advance()
    return g.finish(), resample_uniform(np.array([0.0] + list(T)), np.array(rows).T, step=step)


@pytest.mark.parametrize('end', ['at tend', 'a rounding below a grid point', 'short of tend', 'past tend'])
@pytest.mark.parametrize('block', [1, 200, 1 << 30])
def test_a_hand_made_run_s_edges(end, block, monkeypatch):
    rng = np.random.default_rng(7)
    step, tend = 0.1, 1.0
    last = {'at tend': 1.0, 'a rounding below a grid point': 1.0 - 1e-14,
            'short of tend': 0.83, 'past tend': 1.1}[end]
    T = np.sort(rng.uniform(0.0, last, 40))
    T = np.append(T[T < last], last)
    got, want = _hand_made(T, step, tend, block, monkeypatch)
    assert _bytes(*got) == _bytes(*want)
    assert len(got[0]) == {'at tend': 11, 'a rounding below a grid point': 11, 'short of tend': 9,
                           'past tend': 12}[end]


@pytest.mark.parametrize('seed', range(6))
def test_random_runs_any_block(seed, monkeypatch):
    """Steps from a hundredth of the grid's to ten times it, any block."""
    rng = np.random.default_rng(seed)
    h = 10.0 ** rng.uniform(-2, 1, 300)
    T = np.cumsum(h) * 1e-9
    step = float(T[-1] / rng.integers(20, 400))
    for block in (1, 24 * int(rng.integers(1, 30)), 1 << 30):
        got, want = _hand_made(T, step, float(T[-1]), block, monkeypatch, n=int(rng.integers(1, 6)), seed=seed)
        assert _bytes(*got) == _bytes(*want), block


def test_a_run_too_short_is_refused_as_the_whole_run_is(monkeypatch):
    with pytest.raises(ValueError, match='need at least 3 points to interpolate quadratically, got 2'):
        _hand_made([0.5], 0.1, 0.5, 1, monkeypatch)
