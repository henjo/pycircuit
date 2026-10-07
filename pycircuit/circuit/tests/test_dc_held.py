"""Held nodes in the operating point -- SPICE's `.ic` for a transient
(`DC(pin=...)`; `Transient(ic=...)` without `uic`, in
`test_initial_conditions`) and `.nodeset` (`DC(nodeset=...)`): the SPICE
benchmark plan's stage 3.  A held node is not an unknown of the solve: it
is its value exactly, the rest solves around it, and the hold scales with
the source-stepping factor as SPICE scales a pinned row."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import circuit, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.dcanalysis import DC, _Held
from pycircuit.circuit.elements import VS, R
from pycircuit.circuit.elements_hdl import MosLevel1Hdl


def _divider():
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['v'] = VS('top', gnd, v=10.0)
    c['r1'] = R('top', 'mid', r=1e3)
    c['r2'] = R('mid', gnd, r=2e3)
    return c


def test_a_held_node_is_its_value_and_the_rest_solves_around_it():
    """`mid` held at 3 V: the source still sets `top`, and supplies the
    current r1 carries into the hold; unheld, the divider's own answer."""
    res = DC(_divider(), pin={'mid': 3.0}).solve()
    assert res.v('mid') == 3.0 and res.v('top') == 10.0
    assert res.i('v.plus') == pytest.approx(-(10.0 - 3.0) / 1e3, rel=1e-12, abs=0.0)
    assert DC(_divider()).solve().v('mid') == pytest.approx(20.0 / 3.0, rel=1e-12, abs=0.0)


def test_the_hold_scales_with_the_source_stepping_factor():
    """SPICE's pinned row reads `x_k = srcFact * ic`: a held node enters a
    source-stepping rung at `lam` times its value; the reference node at 0."""
    red = _Held(4, 0, {2: 3.0})
    red.lam = 0.1
    x = red.full(np.array([1.0, 2.0]))
    assert x.tolist() == [0.0, 1.0, 0.1 * 3.0, 2.0]
    assert red.reduce(x).tolist() == [1.0, 2.0]
    F, J = red.system(np.arange(4.0), np.arange(16.0).reshape(4, 4))
    assert F.tolist() == [1.0, 3.0] and J.tolist() == [[5.0, 7.0], [13.0, 15.0]]


@pytest.mark.parametrize('pin, says', [
    ({'top': 5.0}, 'a voltage source or an inductor holds'),
    ({'nowhere': 1.0}, 'not in the circuit'),
    ({gnd: 1.0}, 'reference node'),
], ids=['source', 'unknown', 'reference'])
def test_what_cannot_be_held_is_refused(pin, says):
    with pytest.raises(ValueError, match=says):
        DC(_divider(), pin=pin).solve()


def test_unholdable_names_the_nodes_a_source_holds():
    """What the refusal above says, asked before a solve (the importer
    leaves such an `.ic` out, as SPICE lets the source win)."""
    dc = DC(_divider())
    assert dc.unholdable({'top': 5.0, 'mid': 3.0}) == ['top']
    assert dc.unholdable({'mid': 3.0}) == [] and dc.unholdable(None) == []


def test_pcnr_with_held_nodes_is_refused():
    with pytest.raises(ValueError, match='with held nodes'):
        DC(_divider(), pcnr=True, nodeset={'mid': 1.0}).solve()


def _latch():
    """Two resistor-loaded NMOS inverters, cross-coupled: bistable, and from
    zeros its symmetric metastable point."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=5.0)
    c['r1'] = R('vdd', 'q', r=10e3)
    c['r2'] = R('vdd', 'qb', r=10e3)
    c['m1'] = MosLevel1Hdl('q', 'qb', gnd, gnd, vto=1.0, kp=2e-5, w=10e-6, l=1e-6)
    c['m2'] = MosLevel1Hdl('qb', 'q', gnd, gnd, vto=1.0, kp=2e-5, w=10e-6, l=1e-6)
    return c


def test_a_nodeset_picks_each_state_of_a_latch():
    """From zeros the solve finds the metastable point; a nodeset picks a
    state -- a solve with the nodes held, then one without them -- and the
    answer is an operating point of the unheld circuit (an ordinary solve
    from it stays)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        meta = DC(_latch()).solve()
        assert meta.v('q') == pytest.approx(meta.v('qb'), rel=1e-9)
        for hi, lo in (('q', 'qb'), ('qb', 'q')):
            res = DC(_latch(), nodeset={hi: 5.0, lo: 0.0}).solve()
            assert res.v(hi) > 4.0 and res.v(lo) < 1.0, (hi, res.v(hi), res.v(lo))
            again = DC(_latch()).solve(x0=np.asarray(res.x))
            assert np.allclose(np.asarray(again.x), np.asarray(res.x), rtol=1e-9, atol=1e-12)
