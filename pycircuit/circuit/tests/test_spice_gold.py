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
