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


def test_construction_does_not_rebuild_the_map_per_element():
    """Gate 2+.5-3, counted rather than timed.

    This was a timing test asserting a fitted exponent below 1.5.  It passed
    alone and FAILED inside the full suite on a box at load -- a flake, and a
    flaky test is worse than none because the next person learns to re-run it
    rather than read it.  The mechanism 2+.5 fixed is countable, so count it.

    `add_instance` used to call `update_node_map()`, which rebuilds EVERY
    element's entry: N insertions did O(N^2) element-work before any analysis
    ran.  It now marks the map dirty and one rebuild happens on first read, so
    the rebuild count must not grow with the circuit.
    """
    from pycircuit.circuit import circuit as circuit_mod

    calls = {'n': 0}
    real = circuit_mod.SubCircuit.update_node_map

    def counting(self, *a, **k):
        calls['n'] += 1
        return real(self, *a, **k)

    circuit_mod.SubCircuit.update_node_map = counting
    try:
        counts = []
        for N in (10, 40, 160):
            calls['n'] = 0
            cir = _ladder(N)
            cir.elementnodemap          # force the one lazy rebuild
            counts.append(calls['n'])
    finally:
        circuit_mod.SubCircuit.update_node_map = real

    assert counts[0] == counts[-1], \
        'rebuilds grow with circuit size: %r for N = 10, 40, 160' % counts
    assert counts[-1] <= 2, \
        'expected one rebuild on first read, got %d' % counts[-1]


def test_a_large_circuit_can_be_built_at_all():
    """The point of the item: 800 sections used to take ~100 s."""
    t0 = time.perf_counter()
    cir = _ladder(800)
    elapsed = time.perf_counter() - t0
    assert len(cir.elements) == 1600
    assert elapsed < 15.0, 'building 1600 elements took %.1f s' % elapsed
