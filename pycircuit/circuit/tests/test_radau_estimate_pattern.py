"""Radau's error estimate on C's pattern (`_tran_radau.PATTERN_FROM`,
`_est_pattern`, `_c_times`; loop.c's `est_pat` and `tc_cmul`; roadmap 267).
From `PATTERN_FROM` unknowns, where the run's core holds its pattern, the
estimate's and the early check's ``C F1`` add each row's pattern products in
the pattern's order to +0.0 -- `np.bincount`'s sequential sums -- where
numpy's dense matmul read every n x n entry (71 % of the C loop at 4 000
stages); the C loops make the same sums.  Under it, numpy's matmul as
before.  The values move at rounding level only: a zero entry adds nothing,
the nonzero ones are summed in another order."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_radau,
    _tran_radau_c,
    _tran_radau_tc,
    _tran_stageloop_c,
    _tran_tcloop_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit._tran_radau import _c_times
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.transient import Transient


def _pattern(n, rows, cols):
    pat = np.unique(np.asarray(rows, dtype=np.int64) * n + np.asarray(cols, dtype=np.int64))
    r, c = np.divmod(pat, n)
    return pat, r, c


def test_the_product_is_each_rows_products_summed_in_order():
    """`_c_times` on a pattern: every row's products ``C[i, j] F[j]`` over
    its pattern entries, column by column from +0.0, each addition rounded --
    rows of up to 40 entries (where a pairwise sum would part from it),
    entries of every scale with cancellation, zeros on the pattern, signed
    zeros; without a pattern, numpy's matmul."""
    rng = np.random.default_rng(267)
    for _ in range(40):
        n = int(rng.integers(5, 60))
        Cm = np.zeros((n, n))
        rows, cols = [], []
        for i in range(n):
            k = int(rng.integers(0, min(n, 40) + 1))
            js = rng.choice(n, size=k, replace=False)
            Cm[i, js] = rng.standard_normal(k) * 10.0 ** rng.integers(-15, 3, size=k)
            Cm[i, js[: k // 4]] = 0.0
            rows += [i] * k
            cols += list(js)
        Cm[Cm == 0.0] = rng.choice((0.0, -0.0), size=int(np.count_nonzero(Cm == 0.0)))
        F = rng.standard_normal(n) * 10.0 ** rng.integers(-3, 9, size=n)
        ep = _pattern(n, rows, cols)
        want = np.zeros(n)
        for p, i, j in zip(*ep):
            want[i] = want[i] + Cm[i, j] * F[j]
        got = _c_times(Cm, F, ep)
        assert np.array_equal(got, want) and not np.signbit(got[got == 0.0]).any()
        assert np.array_equal(_c_times(Cm, F, None), Cm @ F)


def test_one_entry_a_row_is_numpys_matmul_bit_for_bit():
    """A row of one nonzero is its product, as numpy's matmul makes it: a
    capacitor to ground a node, the pattern's other entries zero."""
    rng = np.random.default_rng(3)
    n = 300
    Cm = np.diag(rng.random(n) * 1e-14)
    F = rng.standard_normal(n) * 1e9
    rows, cols = np.nonzero(np.ones((n, n)) * (np.abs(np.subtract.outer(np.arange(n), np.arange(n))) <= 2))
    assert np.array_equal(_c_times(Cm, F, _pattern(n, rows, cols)), Cm @ F)


def _chain(n, diodes=True, victim=False):
    """A pulsed MOS chain, a gate-drain capacitor at every stage and one
    three stages back at every third (C's rows 3-5 entries), a junction's
    depletion charge at every fifth node (C moving with the state);
    `victim`: a node coupled to every other stage by 1 fF, 0.1 fF and
    100 kOhm to ground -- the estimate's largest entry on a row of ~20."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-8, tr=1e-8, tf=1e-8, pw=2.4e-7, per=5e-7)
    for k in range(n):
        g = 'g0' if k == 0 else f'd{k - 1}'
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'cm{k}'] = C(f'd{k}', g, c=2e-15)
        if k >= 3 and k % 3 == 0:
            c[f'cl{k}'] = C(f'd{k}', f'd{k - 3}', c=1e-15)
        if diodes and k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    if victim:
        c.add_node('v')
        c['rv'] = R('v', gnd, r=1e5)
        c['cv'] = C('v', gnd, c=1e-16)
        for k in range(0, n, 2):
            c[f'cv{k}'] = C('v', f'd{k}', c=1e-15)
    return c


def _run(monkeypatch, loops, frm, n=40, transform=True, packed='size', tend=3e-7, victim=False):
    """The run's bytes and statistics, its estimate's and early check's
    counts, and the `est_pat` each C binding took."""
    for mod in (_tran_radau_tc, _tran_radau_c, _tran_tcloop_c, _tran_loop_c):
        monkeypatch.setattr(mod, 'ENABLED', loops)
    monkeypatch.setattr(_tran_tcloop_c, 'PACKED', packed)
    monkeypatch.setattr(_tran_radau, 'PATTERN_FROM', frm)
    flags = []
    bind = _tran_stageloop_c._bind_pattern

    def spy(run, core, keep):
        bind(run, core, keep)
        flags.append(int(run.s.est_pat))
    monkeypatch.setattr(_tran_stageloop_c, '_bind_pattern', spy)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(n, victim=victim), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                       radau_transform=transform)
        res = tr.solve(tend=tend, timestep=1e-8)
    got = _paths.since(before)
    st = tr.statistics
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps)
    counts = {k: v for k, v in got.items() if k.startswith(('radau.est:', 'radau.early:'))}
    loop = sum(got.get(k + ':attempts', 0) for k in ('tcloop_c', 'stageloop_c', 'loop_c'))
    return run, counts, flags, loop, np.asarray(res.x, float)


@pytest.mark.parametrize('case', ['transform', 'stage', 'packed', 'victim'])
def test_the_c_loops_make_the_pythons_run_on_the_pattern(monkeypatch, case):
    """With the pattern from 10 unknowns: the C loop's run -- the
    transform's estimate and early check (`tc_estimate`, `tc_check`), the
    stage loop's fresh estimate (`stage_estimate`, the transform off), the
    packed sets' (250 stages) -- is the Python's, bit for bit, its estimates
    and checks counted alike, every binding with `est_pat` set."""
    kw = {'transform': {}, 'stage': {'transform': False},
          'packed': {'n': 250, 'packed': 'always'}, 'victim': {'victim': True}}[case]
    want, wc, _f, none, _x = _run(monkeypatch, False, 10, **kw)
    got, gc, flags, loop, _x = _run(monkeypatch, True, 10, **kw)
    assert none == 0 and loop > 20 and flags and set(flags) == {1}
    assert got == want and gc == wc
    assert sum(v for k, v in wc.items() if k.startswith('radau.est:')) > 20
    if case != 'stage':
        assert wc.get('radau.early:checks', 0) > 0


def test_under_the_size_the_estimate_is_numpys_and_the_switch_is_live(monkeypatch):
    """At the default `PATTERN_FROM` a 45-unknown run takes numpy's matmul
    (no binding sets `est_pat`), its C loop the Python's run; from 10 its
    values move -- at rounding level, the final state within 1e-9 of the
    tolerance."""
    want, _c, _f, _l, xd = _run(monkeypatch, False, _tran_radau.PATTERN_FROM)
    got, _c, flags, loop, _x = _run(monkeypatch, True, _tran_radau.PATTERN_FROM)
    assert got == want and loop > 20 and flags and set(flags) == {0}
    pat, _c, _f, _l, xp = _run(monkeypatch, True, 10)
    assert pat[0] != want[0]
    if xp.shape == xd.shape and pat[1] == want[1]:
        assert np.max(np.abs(xp - xd) / (1e-4 * np.abs(xd) + 1e-6)) < 1e-6
