"""The fetched benchmark decks against Xyce's gold outputs, scored by
Xyce's own metric (`spiceoutput.xyce_verify`; 1 or below passes) -- the
SPICE benchmark plan's stage 4.  Each test skips where the data is not
fetched (`benchmarks/fetch_spice_suite.py`)."""
import warnings

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import circuit, spice_import
from pycircuit.utilities import spiceoutput


def _data(*paths):
    got = [benchdata.spice_data(p) for p in paths]
    if None in got:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    return got


def _score(imp, gold, columns, **solve):
    """Run `imp`'s transient and score each printed column (name -> (node,
    added offset)) against `gold`: {column: (metric, rising crossings ours,
    gold's)}."""
    tr, kw = imp.transient()
    kw.update(solve)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**kw)
    names, g = spiceoutput.read_prn(gold)
    lower = [n.lower() for n in names]
    out = {}
    for col, (node, add) in columns.items():
        w = res.v(node)
        t, v = np.asarray(w.x[0], dtype=float), np.asarray(w.y, dtype=float) + add
        gc = g[:, lower.index(col.lower())]
        out[col] = spiceoutput.xyce_verify(t, v, g[:, 1], gc)
    return out


def test_the_4049_oscillator_passes_xyce_s_own_metric_against_its_gold():
    """Two CD4049UB inverters (MOS level 1, overlap and junction
    capacitances) in an RC relaxation oscillator, started from zeros
    (NOOP): every printed column within Xyce's tolerance of Xyce's
    waveform (measured 2026-10-07: 0.358, 0.385, 0.366, 0.359 -- the
    rising edges within 0.01 us of the gold's over a ~197 us period)."""
    deck, gold = _data('Netlists/4049OSC/4049osc.cir', 'OutputData/4049OSC/4049osc.cir.prn')
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    got = _score(imp, gold, {'{V(8)+4}': ('8', 4.0), '{v(5)+4}': ('5', 4.0),
                             '{V(1)+4}': ('1', 4.0), '{V(3)+4}': ('3', 4.0)})
    assert all(e < 1.0 for e in got.values()), got
    assert max(got.values()) < 0.5, got


def _rising(t, y, level):
    """Interpolated upward crossings of `level`."""
    k = np.nonzero((y[:-1] < level) & (y[1:] >= level))[0]
    return t[k] + (level - y[k]) * (t[k + 1] - t[k]) / (y[k + 1] - y[k])


def test_the_4049_oscillator_s_periodic_steady_state_is_the_gold_s_settled_cycle():
    """Stage 8: the 4049 oscillator's limit cycle by autonomous shooting --
    seeded as the deck starts (NOOP: zeros, two periods of its own
    transient), on a grid from one adaptive gear period (`lte_grid`), the
    period an unknown seeded at 197 us.  Against the gold's LAST cycle
    (281.74 .. 478.73 us, 196.9905 us), phase-aligned at v(3)'s rising
    edge: measured 2026-10-09 the period 196.98892 us and every node
    within 3.0 mV over swings of 5-10 V (radau's own grid: 196.99028 us,
    2.4 mV; trbdf2's 196.99004 us, 1.0 mV -- a uniform grid of 400 is
    0.27 V off at the switching edge, the derived grid is the point)."""
    from pycircuit.circuit.shooting.pss import PSS
    deck, gold = _data('Netlists/4049OSC/4049osc.cir', 'OutputData/4049OSC/4049osc.cir.prn')
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    cir = imp.circuit
    T0 = 197e-6
    tr, kw = imp.transient()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        start = np.asarray(tr.solve(tend=2 * T0, timestep=kw['timestep']).x)[:, -1]
        pss = PSS(cir, method='gear', epar=imp.epar(), reltol=imp.options['reltol'])
        fracs, seed = pss.lte_grid(T0, x0=start, tstab=2 * T0, fold=False)
        out = pss.solve(period=pss.lte_period, grid=fracs, x0=seed, maxiterations=30)
    assert pss.converged
    names, g = spiceoutput.read_prn(gold)
    tg = g[:, names.index('TIME')]
    g0, g1 = _rising(tg, g[:, names.index('{V(3)+4}')], 6.5)[-2:]
    ## (8e-6 measured; the gold's own cycles still move 3e-4 from the second
    ## to the third, a uniform grid of 400 is 1e-3 off)
    assert abs(pss.period / (g1 - g0) - 1.0) < 2e-5, (pss.period, g1 - g0)
    T = pss.period

    def wave(node):
        w = out['tpss'].v(node)
        t, y = np.asarray(w.x[0], dtype=float), np.asarray(w.y, dtype=float)
        ## (three periods laid end to end, for any phase)
        return np.r_[t, t[1:] + T, t[1:] + 2 * T], np.r_[y, y[1:], y[1:]]

    t3, y3 = wave('3')
    c0 = _rising(t3, y3, 2.5)[0]
    sel = (tg >= g0) & (tg <= g1)
    worst = {}
    for col, node in (('{V(8)+4}', '8'), ('{V(5)+4}', '5'), ('{V(1)+4}', '1'),
                      ('{V(3)+4}', '3')):
        t, y = wave(node)
        d = np.interp(tg[sel] - g0 + c0, t, y) - (g[sel, names.index(col)] - 4.0)
        worst[node] = float(np.max(np.abs(d)))
    assert max(worst.values()) < 0.01, worst


def _driven_pss(deck, period, method, steps=200):
    """A driven deck's periodic steady state seeded at its operating point,
    as harmonic balance starts (the decks' `.op`; from zeros the
    common-emitter stage's shooting Newton diverges -- measured)."""
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.shooting.pss import PSS
    imp = spice_import.import_netlist(deck)
    cir = imp.circuit
    xdc = np.asarray(DC(cir, epar=imp.epar()).solve().x, dtype=float).ravel()
    pss = PSS(cir, method=method, epar=imp.epar(), **imp.options)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = pss.solve(period=period, timestep=period / steps, maxiterations=40,
                        x0=np.delete(xdc, cir.get_node_index(circuit.gnd)))
    assert pss.converged
    return out


def _hb_scores(out, td, fd, columns, period, harmonics):
    """Each column (name -> (node, node or None)) against Xyce's HB gold:
    {column: (xyce_verify on the time points, worst harmonic error over
    `harmonics` relative to the fundamental)}.  Xyce's frequency output is
    the two-sided Fourier coefficient; `fpss` is the RMS phasor, sqrt(2)
    times it above DC (on the 1 V source: -0.7071j against -0.5j)."""
    names, g = spiceoutput.read_prn(td)
    names = [n.upper() for n in names]
    fnames, gf = spiceoutput.read_prn(fd)
    fnames = [n.upper() for n in fnames]
    got = {}
    for col, (a, b) in columns.items():
        w = out['tpss'].v(a, b) if b else out['tpss'].v(a)
        ours = np.interp(g[:, 1], np.asarray(w.x[0], dtype=float), np.asarray(w.y, dtype=float))
        metric = spiceoutput.xyce_verify(g[:, 1], ours, g[:, 1], g[:, names.index(col)])
        X = gf[:, fnames.index(f'RE({col})')] + 1j * gf[:, fnames.index(f'IM({col})')]
        w = out['fpss'].v(a, b) if b else out['fpss'].v(a)
        fx, fy = np.asarray(w.x[0], dtype=float), np.asarray(w.y)
        err = []
        for h in harmonics:
            mine = fy[np.argmin(np.abs(fx - h / period))] / (1.0 if h == 0 else np.sqrt(2))
            err.append(abs(mine - X[np.argmin(np.abs(gf[:, 1] - h / period))]))
        x1 = abs(X[np.argmin(np.abs(gf[:, 1] - 1.0 / period))])
        got[col] = (metric, max(err) / x1)
    return got


def test_the_common_emitter_stage_s_periodic_steady_state_is_xyce_s_harmonic_balance():
    """Stage 8, driven: a 2N2222 common-emitter stage driven at 1 MHz into
    clipping, by radau shooting on 200 steps, against Xyce's harmonic
    balance (50 harmonics).  Measured 2026-10-09: xyce_verify 0.094 (ve)
    and 0.808 (out); harmonics 0..5 within 2e-4 of the fundamental (h1 of
    out: 0.15509+2.6024j against 0.15505+2.6023j).  The time-domain
    residual of `out` (12 mV over a 10 V swing at 1600 steps, at the
    clipping edge) is the GOLD's: its harmonics 49-50 hold 1.9 mV where
    ours hold 0.7, and our orbit carries 3.4 mV rms above the 50th --
    what a 50-harmonic solve cannot represent."""
    R = 'Netlists/HB/common_emitter_hb.cir'
    deck, td, fd = _data(R, 'OutputData/HB/common_emitter_hb.cir.HB.TD.prn',
                         'OutputData/HB/common_emitter_hb.cir.HB.FD.prn')
    circuit.default_toolkit = circuit.numeric
    out = _driven_pss(deck, 1e-6, 'radau')
    got = _hb_scores(out, td, fd, {'V(VE)': ('ve', None), 'V(OUT)': ('out', None)},
                     1e-6, range(6))
    assert all(m < 1.0 for m, _e in got.values()), got
    assert got['V(VE)'][0] < 0.2, got
    assert all(e < 1e-3 for _m, e in got.values()), got


def test_the_gilbert_cell_s_periodic_steady_state_is_xyce_s_harmonic_balance():
    """Stage 8, driven: a bipolar Gilbert cell (six QB2T2222 with irb, PTF
    and a substrate junction) with a 10 kHz input, by gear shooting on 200
    steps, against Xyce's harmonic balance (20 harmonics).  Measured
    2026-10-09: the output V(5,3) within 0.17 mV over 2.47 V (xyce_verify
    0.016), the input 0.009; h1 -0.63630j against -0.63632j, h3 and h5 to
    4 digits, the even harmonics zero in both (a balanced cell).  Gear in
    0.5 s; radau the same answer in 1.2 s (222 s until `PSS.solve` held
    BLAS at one thread, 2026-10-09)."""
    R = 'Netlists/HB/gilbert_cell_hb.cir'
    deck, td, fd = _data(R, 'OutputData/HB/gilbert_cell_hb.cir.HB.TD.prn',
                         'OutputData/HB/gilbert_cell_hb.cir.HB.FD.prn')
    circuit.default_toolkit = circuit.numeric
    out = _driven_pss(deck, 1e-4, 'gear')
    got = _hb_scores(out, td, fd, {'V(15,10)': ('15', '10'), 'V(5,3)': ('5', '3')},
                     1e-4, range(6))
    assert all(m < 0.1 for m, _e in got.values()), got
    assert all(e < 1e-3 for _m, e in got.values()), got


def test_the_ring_oscillator_s_periodic_steady_state():
    """Stage 9: the CircuitSim90 ring (MOS level 2, nine stages) by
    autonomous shooting -- settled by 30 ns of its own transient from the
    deck's `.ic`, a grid from one adaptive period on the GLOBAL reference
    (`lte_grid`'s per-unknown default drew 24k points through the fast
    edges), gear, the supply currents out of the error test (as the suite
    runs this deck).  Measured 2026-10-10: the period converges with the
    grid, 10.877 / 10.856 / 10.852 ns at 320 / 1013 / 2864 points; the
    gold's cycles are 10.56 and ngspice's on our deck 10.68 -- the 1.6 %
    to ngspice is the gate charge's open gap (doc/hdl_roadmap_260824.md,
    84), bounded here so it cannot grow unseen."""
    from pycircuit.circuit.shooting.pss import PSS
    deck, = _data('Netlists/CircuitSim90/MOS2/ring.cir')
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    tr, kw = imp.transient(lte_algebraic='skip')
    T0 = 10.56e-9
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        start = np.asarray(tr.solve(tend=30e-9, timestep=kw['timestep']).x)[:, -1]
        pss = PSS(imp.circuit, method='gear', epar=imp.epar(), lte_algebraic='skip',
                  **imp.options)
        fracs, seed = pss.lte_grid(T0, x0=start, tstab=2 * T0, fold=False,
                                   relref='sigglobal', reltol=1e-4)
        pss.solve(period=pss.lte_period, grid=fracs, x0=seed, maxiterations=40)
    assert pss.converged
    assert abs(pss.period / 10.852e-9 - 1.0) < 5e-3, pss.period     # its grid-converged value
    assert abs(pss.period / 10.68e-9 - 1.0) < 0.025, pss.period     # ngspice: the open gap


def test_gm17_s_supply_currents_in_the_error_test_force_few_steps():
    """The step controller's retry order (2026-10-10): gm17 (CircuitSim90,
    MOS level 3) with its supply currents JUDGED -- 'judge', the error
    test's default until 2026-10-13, not the suite's skip then.  Their
    estimates fall as ~h^1.1-1.4
    between retries at a point; the cube-law retries undershot and 26
    steps were force-accepted over tolerance (and until the gate charge's
    onset split, 10700247, the run stopped at 5.75 ns).  With the second
    retry at the order its tries measured: 2."""
    deck, = _data('Netlists/CircuitSim90/MOS3/gm17.cir')
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(deck)
    tr, kw = imp.transient(lte_algebraic='judge')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(**kw)
    assert tr.statistics.force_accepts <= 5, tr.statistics.force_accepts


@pytest.mark.parametrize('radau', [False, True])
def test_gm3_s_supply_current_by_the_charge_it_moves_forces_no_step(radau):
    """`lte_algebraic='charge'`, the default since 2026-10-13
    (doc/hdl_roadmap_260824.md, 76, 129 and 130): gm3 (CircuitSim90, MOS
    level 3) with its supply's current judged as every unknown force-accepts
    steps over tolerance -- the current is the gate charges' slope (it closes
    a loop of capacitors), and its estimate does not fall with the step;
    judged by the charge it moves, the multistep error test forces none, in
    fewer steps (measured 854 / 7 forced -> 631 / 0).  Radau gains too (593
    / 40 -> 435 / 23); its own step controller -- an RMS norm, fixed-order
    retries -- still forces some at the corners."""
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    deck, = _data('Netlists/CircuitSim90/MOS3/gm3.cir')
    circuit.default_toolkit = circuit.numeric
    got = {}
    for mode in ('judge', 'charge'):
        imp = spice_import.import_netlist(deck)
        extra = {'integrator': RadauIIA3Integrator()} if radau else {}
        tr, kw = imp.transient(lte_algebraic=mode, **extra)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr.solve(**kw)
        got[mode] = (tr.statistics.accepted_steps, tr.statistics.force_accepts)
    assert got['judge'][1] > 0 and got['charge'][0] < got['judge'][0], got
    if radau:
        assert got['charge'][1] < got['judge'][1], got
    else:
        assert got['charge'][1] == 0, got
