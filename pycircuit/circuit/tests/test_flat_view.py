"""A hierarchical circuit is evaluated through its flat view
(`SubCircuit._flat_view`): its leaf elements on one level, its own unknowns
in its own order -- so the C paths, which take a flat list of leaves and
declined a subcircuit instance, serve it; its results report against it."""
import warnings

import numpy as np
import pytest
from sympy import Symbol

from pycircuit.circuit import _paths, circuit, gnd, symbolic
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.dcanalysis import DC, DCSweep
from pycircuit.circuit.elements import VS, C, R
from pycircuit.circuit.transient import Transient
from pycircuit.utilities.param import Parameter

TK = circuit.numeric


class _Stage(SubCircuit):
    terminals = ('a', 'b')
    instparams = [Parameter('rload', default=1e3)]  # noqa: RUF012

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self['R1'] = R('a', 'm', r='rload', toolkit=TK)
        self['D1'] = eh.DiodeSpiceHdl('m', 'b', IS=1e-14, toolkit=TK)
        self['C1'] = C('m', 'b', c=1e-9, toolkit=TK)


def _hier():
    c = SubCircuit(toolkit=TK)
    n1, n2 = c.add_nodes('in', 'mid')
    c['V'] = VS(n1, gnd, v=2.0, toolkit=TK)
    c['X1'] = _Stage(n1, n2, toolkit=TK)
    c['X2'] = _Stage(n2, gnd, rload=2e3, toolkit=TK)
    return c


def _run(c, counts=None):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        b = _paths.snapshot()
        x = np.asarray(DC(c).solve().x, dtype=float)
        r = Transient(c, outputstep=1e-7).solve(tend=2e-6, timestep=1e-8, x0=x)
        if counts is not None:
            counts.update(_paths.since(b))
    return x, np.asarray(r.v('mid', gnd).get_y(), dtype=float)


def test_a_hierarchical_circuit_is_served_by_the_c_paths(monkeypatch):
    """Through the view the PCNR C solve and the C Newton serve a circuit of
    subcircuit instances, and its answers are the hierarchical evaluation's
    (until 2026-10-22 both declined it -- 'unservable' -- and every element
    ran the Python paths, 7..40x slower)."""
    got = {}
    x, y = _run(_hier(), got)
    assert got.get('pcnr_c:served', 0) > 0 and got.get('newton_c:served', 0) > 0
    assert not any(k.endswith(':unservable') for k, v in got.items() if v)
    monkeypatch.setattr(SubCircuit, '_flat_view', lambda self: self)
    old = {}
    x0, y0 = _run(_hier(), old)
    ## (without the view PCNR does not see the diodes at all -- they sit in
    ## the instances -- and the C Newton declines: nothing served by C)
    assert old.get('pcnr_c:served', 0) == 0 and old.get('newton_c:served', 0) == 0
    assert np.allclose(x, x0, rtol=1e-8, atol=1e-12)
    assert len(y) == len(y0) and np.allclose(y, y0, rtol=1e-8, atol=1e-10)


def test_the_view_keeps_the_circuit_s_unknowns_and_its_results_report_against_it():
    c = _hier()
    dc = DC(c)
    assert dc.circuit is c and dc.cir is not c
    assert dc.cir.nodes == c.nodes and dc.cir.branches == c.branches
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = dc.solve()
    assert res.circuit is c
    assert np.isfinite(float(res.v('X1.m', gnd)))
    assert dc.cir['X1.R1'] is c['X1']['R1']


def test_a_subcircuit_s_parameter_swept_reaches_the_view():
    """DCSweep writes the instance's own parameter into the circuit given;
    the view's leaves resolve it there."""
    vals = [1e3, 2e3, 4e3]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = DCSweep(_hier()).solve('X2', 'rload', vals)
        swept = np.asarray(res.v('mid', gnd), dtype=float)
        single = []
        for v in vals:
            c = _hier()
            c['X2'].ipar.rload = v
            single.append(float(DC(c).solve().v('mid', gnd)))
    ## (each point converged to the DC's own tolerance, reltol 1e-4, from
    ## the last point or from zeros)
    assert np.allclose(swept, single, rtol=2e-4, atol=1e-9)
    assert len(set(np.round(swept, 6))) == 3


def test_a_topology_change_after_the_analysis_was_made_is_in_its_solve():
    c = _hier()
    dc = DC(c)
    c['X1']['R2'] = R('m', 'b', r=1e4, toolkit=TK)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = np.asarray(dc.solve().x, dtype=float)
        want = np.asarray(DC(c).solve().x, dtype=float)
    assert got.tobytes() == want.tobytes()


def test_a_leaf_keeps_its_own_toolkit():
    """A circuit may mix toolkits -- symbolic devices in a numeric circuit
    (the regulated-cascode example): building the view leaves each leaf's
    toolkit as it was."""
    from pycircuit.circuit.mos import MOS
    c = SubCircuit(toolkit=TK)
    nin, nout = c.add_nodes('in', 'out')
    c['M'] = MOS(nin, nout, gnd, gnd, gm=Symbol('gm'), toolkit=symbolic)
    leaf = c['M']['Igm']
    assert leaf.toolkit is symbolic
    DC(c)
    assert leaf.toolkit is symbolic


@pytest.mark.parametrize('nested', [False, True], ids=['one level', 'two levels'])
def test_a_flat_circuit_gets_no_view_and_a_nested_one_a_complete_one(nested):
    flat = SubCircuit(toolkit=TK)
    flat.add_node('a')
    flat['V'] = VS('a', gnd, v=1.0, toolkit=TK)
    flat['R'] = R('a', gnd, r=1e3, toolkit=TK)
    assert flat._flat_view() is flat
    c = _hier()
    if nested:
        outer = SubCircuit(toolkit=TK)
        outer.add_nodes('in', 'mid')
        outer['V'] = VS('in', gnd, v=1.0, toolkit=TK)
        top = _Stage('in', 'mid', toolkit=TK)
        outer['T'] = top
        outer['Y'] = _Stage('mid', gnd, toolkit=TK)
        c = outer
    v = c._flat_view()
    assert v is not c and not any(isinstance(e, SubCircuit) for e in v.elements.values())
    assert v.nodes == c.nodes and v.branches == c.branches
