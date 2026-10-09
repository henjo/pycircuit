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
