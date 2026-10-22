"""The dense multistep walk carries its monodromy on the live columns only
(`_steps._LiveMonodromy`): a column zero in the whole carried state stays
zero, so it is dropped and put back at the end; the columns left are those
a charge reaches.  Each column is computed as it was."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import C, Diode, R, VSin
from pycircuit.circuit.linearsolver import DenseSolver
from pycircuit.circuit.shooting import PSS, _steps
from pycircuit.circuit.shooting._steps import _LiveMonodromy, _lmm_recursion
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._shooting_fixtures import _q20_rlc
from pycircuit.circuit.tests._warnpolicy import quiet

MAPS = [('trap', 'plain', 1), ('euler', 'plain', 1), ('gear2', 'pair', 2)]


def _clipper():
    """A diode clipper behind a resistive divider: two capacitive nodes of
    five unknowns (the source's branch and three algebraic nodes beside)."""
    c = SubCircuit()
    c.add_nodes('a', 'b', 'd')
    c['vs'] = VSin('a', gnd, va=2.0, freq=1e3)
    c['R1'] = R('a', 'b', r=1e3)
    c['D'] = Diode('b', gnd)
    c['C1'] = C('b', gnd, c=1e-7)
    c['R2'] = R('b', 'd', r=2e3)
    c['R3'] = R('d', gnd, r=2e3)
    c['C2'] = C('d', gnd, c=1e-8)
    return c


class _Counting(DenseSolver):
    def __init__(self):
        self.rhs = []

    def solve(self, A, b, toolkit):
        self.rhs.append(np.shape(b))
        return DenseSolver.solve(self, A, b, toolkit)


def _walk(make, method, kind, width, live, monkeypatch, solver=None):
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_steps, 'LIVE_COLUMNS', live)
    cir = make()
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-9, linearsolver=solver)
    pss._open_at_x0 = False
    pss.autonomous = False
    assert pss._map_kind() == kind
    times, hs = pss._period_grid(1e-3, 60, None)
    with quiet():
        return pss._walk(kind, np.full(width * m, 0.1), times, hs,
                         T=1e-3).monodromy(), m


@pytest.mark.parametrize('make', [_q20_rlc, _clipper], ids=['rlc', 'clipper'])
@pytest.mark.parametrize(('method', 'kind', 'width'), MAPS, ids=[m[0] for m in MAPS])
def test_the_monodromy_on_the_live_columns_is_every_column_s(make, method, kind, width, monkeypatch):
    full, m = _walk(make, method, kind, width, False, monkeypatch)
    solver = _Counting()
    got, _ = _walk(make, method, kind, width, True, monkeypatch, solver)
    assert got.shape == full.shape == (width * m, width * m)
    assert np.array_equal(got != 0, full != 0)
    assert np.allclose(got, full, rtol=1e-13, atol=0)
    ## the steps after the opening solve the live columns alone: those of
    ## the nodes and branches a charge reaches (two of five here)
    widths = [r[1] for r in solver.rhs if len(r) == 2]
    assert widths[0] == width * m and widths[-1] == 2 * width


@pytest.mark.parametrize('make', [_q20_rlc, _clipper], ids=['rlc', 'clipper'])
@pytest.mark.parametrize('method', ['radau', 'trbdf2', 'esdirk43'])
def test_a_stage_walk_s_monodromy_on_the_live_columns_is_every_column_s(make, method, monkeypatch):
    """`_LiveColumns`: the stage family's walk, through the C step map where
    it serves (`_sens_c`; the RLC's Radau) and the Python step elsewhere
    (a junction: the clipper's diode)."""
    full, m = _walk(make, method, 'stage', 1, False, monkeypatch)
    before = _paths.snapshot()
    got, _ = _walk(make, method, 'stage', 1, True, monkeypatch)
    served = _paths.since(before).get('pss.sens:served', 0)
    assert got.shape == full.shape == (m, m)
    assert np.array_equal(got != 0, full != 0)
    assert np.allclose(got, full, rtol=1e-13, atol=0)
    if make is _q20_rlc and method == 'radau':
        assert served > 0


@pytest.mark.parametrize('method', ['trap', 'gear2', 'radau'])
def test_a_pss_on_the_live_columns_answers_as_on_every_column(method, monkeypatch):
    circuit.default_toolkit = circuit.numeric
    out = {}
    for live in (False, True):
        monkeypatch.setattr(_steps, 'LIVE_COLUMNS', live)
        ## (the orbit's resolution is not the question: both runs read the
        ## same grid)
        with quiet(), warnings.catch_warnings():
            warnings.simplefilter('ignore', AccuracyWarning)
            res = PSS(_clipper(), method=method).solve(period=1e-3, timestep=2e-6)
        out[live] = np.asarray(res['tpss'].v('d').get_y(), dtype=float)
    assert np.allclose(out[True], out[False], rtol=1e-12, atol=1e-15)


def test_a_non_finite_step_runs_at_full_width():
    """A dropped column multiplied into a non-finite `C` is NaN at full
    width (``0 * inf``), not zero: that step is recomputed from the full
    state, the NaN where every column put it, and the walk carries every
    column from there on."""
    rng = np.random.default_rng(3)
    m = 4
    Jf = np.eye(m) * 3 + rng.standard_normal((m, m)) * 0.1
    Cfin = np.zeros((m, m))
    Cfin[0, 0], Cfin[1, 1], Cfin[0, 1] = 1.0, 2.0, 0.5

    def solve(S):
        return np.linalg.solve(Jf, S)
    alphas, b = (1.5, -2.0, 0.5), 0.0
    eye = np.eye(m)
    mono = _LiveMonodromy(eye, eye, np.zeros((m, m)))
    ref = [eye, eye]
    refq = np.zeros((m, m))
    Cs = [Cfin, Cfin]
    for _ in range(3):
        mono.step(Cs, Cfin, alphas, b, solve)
        p, refq = _lmm_recursion(ref, Cs, refq, Cfin, alphas, b, solve)
        ref = [p, ref[0]]
    assert mono.cols is not None and len(mono.cols) == 2
    Cbad = Cfin.copy()
    Cbad[2, 3] = np.inf
    before = _paths.snapshot()
    with np.errstate(invalid='ignore'):
        mono.step([Cbad, Cfin], Cfin, alphas, b, solve)
        p, refq = _lmm_recursion(ref, [Cbad, Cfin], refq, Cfin, alphas, b, solve)
    ref = [p, ref[0]]
    assert _paths.since(before).get('pss.live:nonfinite') == 1
    assert mono.off and mono.cols is None
    got = mono.ring()
    assert np.isnan(ref[0]).any()
    ## (the NaN where every column puts it; the finite entries the full
    ## width's to the rounding -- the products read `C`'s nonzero columns)
    for g, r in zip(got, ref):
        assert np.array_equal(np.isnan(g), np.isnan(r))
        assert np.allclose(np.nan_to_num(g), np.nan_to_num(r), rtol=1e-13, atol=0)
