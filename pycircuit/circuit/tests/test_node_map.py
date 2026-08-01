"""Tests for the lazily-rebuilt node map (stage 2+.5)."""
import time

import numpy as np
import pytest

from pycircuit.circuit import numeric, gnd, SubCircuit
from pycircuit.circuit.elements import R, C, L, VS, IS, Diode


def _ladder(N):
    cir = SubCircuit(toolkit=numeric)
    for i in range(N):
        cir.add_node('n%d' % i)
    cir['vs'] = VS('n0', gnd, v=1.0)
    for i in range(N - 1):
        cir['R%d' % i] = R('n%d' % i, 'n%d' % (i + 1), r=100.0 * (i + 1))
    for i in range(N):
        cir['C%d' % i] = C('n%d' % i, gnd, c=1e-9 * (i + 1))
    return cir


def test_node_indices_are_unchanged_by_the_lazy_rebuild():
    """The map is the MATRIX ROW INDEX; reordering it permutes every G, C and J.

    Checked against indices derived independently from `cir.nodes` rather than
    against a recorded snapshot, so the test states the INVARIANT (an element's
    entry is the position of its nodes in `cir.nodes`, followed by its branch
    rows) instead of pinning whatever the code produced on the day.
    """
    cir = _ladder(12)
    cir.update_iparv()
    n_nodes = len(cir.nodes)
    for instance, element in cir.elements.items():
        nodemap = list(map(int, cir.elementnodemap[instance]))
        term_map = cir.term_node_map[instance]
        expected = [cir.nodes.index(term_map[t]) for t in element.terminals]
        expected += [cir.nodes.index(nd)
                     for nd in element.non_terminal_nodes(instance)]
        assert nodemap[:len(expected)] == expected, \
            '%s: node rows %r != %r' % (instance, nodemap[:len(expected)], expected)
        for branch_row in nodemap[len(expected):]:
            assert branch_row >= n_nodes, \
                '%s: branch row %d overlaps the node rows' % (instance, branch_row)


def test_first_occurrence_wins_like_list_index():
    """`list.index` returns the FIRST match; a dict comprehension keeps the last.

    If two Node objects in `cir.nodes` ever compare equal, the naive rewrite would
    silently renumber matrix rows.  `Node.__eq__` is name-based, so two equal-named
    nodes are exactly the case that would bite.
    """
    from pycircuit.circuit.circuit import Node
    nodes = [Node('a'), Node('b'), Node('a')]
    index = {}
    for i, nd in enumerate(nodes):
        if nd not in index:
            index[nd] = i
    assert index[Node('a')] == 0 == nodes.index(Node('a'))
    naive = {nd: i for i, nd in enumerate(nodes)}
    assert naive[Node('a')] == 2, 'the naive form is what this guards against'


def test_map_is_rebuilt_after_elements_are_added():
    """Deferral must not mean staleness -- every read sees a current map."""
    cir = SubCircuit(toolkit=numeric)
    cir.add_node('a')
    cir['R1'] = R('a', gnd, r=1e3)
    first = {k: list(map(int, v)) for k, v in cir.elementnodemap.items()}
    assert set(first) == {'R1'}

    ## Adding an element introduces a node AND shifts every branch row, which is
    ## why the map cannot be appended to incrementally.
    cir.add_node('b')
    cir['V1'] = VS('b', gnd, v=1.0)
    second = cir.elementnodemap
    assert set(second) == {'R1', 'V1'}
    n_nodes = len(cir.nodes)
    for row in second['V1']:
        assert 0 <= int(row) < cir.n
    assert any(int(row) >= n_nodes for row in second['V1']), \
        'the voltage source should own a branch row'


def test_map_is_rebuilt_after_an_element_is_deleted():
    cir = _ladder(6)
    before = set(cir.elementnodemap)
    del cir['R2']
    after = set(cir.elementnodemap)
    assert 'R2' in before and 'R2' not in after


def test_construction_is_no_longer_quadratic():
    """Gate 2+.5-3: the EXPONENT, not a speedup ratio.

    A constant-factor win would leave the wall in the same place.  Measured
    before this change: N^2.23, 18.9 s for a 400-section ladder and a 1600-section
    one unbuildable inside ten minutes -- which is what stopped stage 7b measuring
    its n=2000 case.

    The bar is deliberately loose (1.5, against a measured 1.15) because this is a
    timing test on a shared box; it is here to catch a return to quadratic, not to
    police a few percent.
    """
    sizes = (50, 100, 200, 400)
    times = []
    for N in sizes:
        t0 = time.perf_counter()
        _ladder(N)
        times.append(max(time.perf_counter() - t0, 1e-6))
    exponent = np.polyfit(np.log(np.array(sizes, dtype=float)),
                          np.log(np.array(times)), 1)[0]
    assert exponent < 1.5, \
        'construction is scaling as N^%.2f; it was N^2.23 before 2+.5' % exponent


def test_a_large_circuit_can_be_built_at_all():
    """The point of the item: 800 sections used to take ~100 s."""
    t0 = time.perf_counter()
    cir = _ladder(800)
    elapsed = time.perf_counter() - t0
    assert len(cir.elements) == 1600
    assert elapsed < 15.0, 'building 1600 elements took %.1f s' % elapsed
