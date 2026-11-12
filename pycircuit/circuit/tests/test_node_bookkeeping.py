"""A circuit's node list (`Circuit.append_node`): membership by name, as
`Node.__eq__` decides it, through a set of the list's names kept beside it --
building a circuit (and the operating point's PCNR view, every run) compares
its nodes O(n) times, where a scan of the list made it O(n^2) -- and the set
dropped wherever the list is edited in place."""
from pycircuit.circuit.circuit import Node, SubCircuit, gnd
from pycircuit.circuit.elements import R


def _ladder(n):
    c = SubCircuit()
    c.add_node('n0')
    for k in range(n):
        c.add_node(f'n{k + 1}')
        c[f'r{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
    c['rg'] = R(f'n{n}', gnd, r=1e3)
    return c


def test_building_a_circuit_compares_its_nodes_linearly(monkeypatch):
    """A ladder of 400 and of 800 resistors: the node comparisons building
    it about double, a few a node -- a scan of the list made them grow with
    the square (80 000 and 320 000)."""
    calls = [0]
    eq = Node.__eq__

    def counted(self, a):
        calls[0] += 1
        return eq(self, a)
    monkeypatch.setattr(Node, '__eq__', counted)
    counts = []
    for n in (400, 800):
        calls[0] = 0
        c = _ladder(n)
        counts.append(calls[0])
        assert len(c.nodes) == n + 2 and len({nd.name for nd in c.nodes}) == n + 2
    assert counts[1] < 2.5 * counts[0] and counts[1] < 20 * 800


def test_a_node_of_a_name_present_is_not_appended_again():
    """Another `Node` of a name the list holds is not appended (`__eq__` is
    by name); the name's entry in `nodenames` becomes it, as before."""
    c = _ladder(3)
    again = Node('n2')
    before = list(c.nodes)
    c.append_node(again)
    assert c.nodes == before and c.nodenames['n2'] is again
    c.append_node(Node('fresh'))
    assert c.nodes[-1] == Node('fresh') and len(c.nodes) == len(before) + 1


class _Two(SubCircuit):
    terminals = ('a', 'b')


def test_an_edit_in_place_drops_the_names():
    """`connect_terminals` replaces a terminal's node -- the list's length
    the same -- and a name it removed is appended again after it."""
    c = _Two()
    c.append_node(Node('x'))
    c.connect_terminals(a=Node('n1'))
    assert [nd.name for nd in c.nodes] == ['b', 'x', 'n1']
    c.append_node(Node('a'))
    assert [nd.name for nd in c.nodes] == ['b', 'x', 'n1', 'a']


def test_a_removal_makes_the_names_again():
    """An element's internal node removed with it (`del`) -- the list
    shorter: its names made again -- the name is appended where it comes
    back."""
    class Pair(SubCircuit):
        terminals = ('p', 'm')

        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.add_node('mid')
            self['r1'] = R('p', 'mid', r=1.0)
            self['r2'] = R('mid', 'm', r=1.0)
    c = SubCircuit()
    c.add_nodes('in', 'out')
    c['x1'] = Pair('in', 'out')
    assert 'x1.mid' in [nd.name for nd in c.nodes]
    del c['x1']
    assert 'x1.mid' not in [nd.name for nd in c.nodes]
    c.append_node(Node('x1.mid'))
    assert [nd.name for nd in c.nodes].count('x1.mid') == 1


def test_a_copy_keeps_its_own_names():
    """A copy's node list is its own: a node appended to the original after
    it is appended to the copy too where it is added there (the names kept
    are the list's they were made from)."""
    import copy
    c = _ladder(3)
    c.append_node(Node('n1'))
    c2 = copy.copy(c)
    assert c2.nodes is not c.nodes
    c.append_node(Node('late'))
    c2.append_node(Node('late'))
    assert [nd.name for nd in c2.nodes].count('late') == 1
    assert [nd.name for nd in c.nodes].count('late') == 1


def test_a_list_put_in_place_is_the_one_its_names_come_from():
    """A node list assigned in place of the circuit's (as `CircuitProxy`
    shares its circuit's) -- of the same length -- is the one membership
    reads: a name it holds is not appended again."""
    c = _ladder(1)
    c.append_node(Node('n1'))
    c.nodes = [Node('q'), Node('r'), Node('s')]
    c.append_node(Node('q'))
    assert [nd.name for nd in c.nodes] == ['q', 'r', 's']


def _scan_index(c, name):
    """`get_node_index`'s answer by the list scan it made before the map."""
    nd = Node(name)
    return c.nodes.index(nd) if nd in c.nodes else None


def _all_agree(c, extra=()):
    """Every node name of `c` (and the names `extra`) found where the scan
    finds it; a name absent raising as the scan's did."""
    import pytest
    for name in [nd.name for nd in c.nodes] + list(extra):
        want = _scan_index(c, name)
        if want is None:
            with pytest.raises(ValueError, match='is not in circuit node list'):
                c.get_node_index(name)
        else:
            assert c.get_node_index(name) == want, name


def test_a_nodes_index_is_found_by_its_name_in_constant_comparisons(monkeypatch):
    """`get_node_index` through the map of the list's places by name kept
    beside it (`_node_index_of`): on a ladder of 4 000 nodes, after the map is
    made, 200 lookups compare two nodes each at most -- the scan compared
    n / 2 on average (~400 000 for these)."""
    c = _ladder(4000)
    c.get_node_index('n1')
    calls = [0]
    eq = Node.__eq__

    def counted(self, a):
        calls[0] += 1
        return eq(self, a)
    monkeypatch.setattr(Node, '__eq__', counted)
    for k in range(0, 4000, 20):
        assert c.get_node_index(f'n{k}') == k
    assert calls[0] <= 2 * 200, calls[0]


def test_the_index_follows_every_edit_of_the_list():
    """The lookups the scan's answers through each edit of the list: a node
    moved (`add_terminals`), replaced in place (`connect_terminals`), removed
    (`del`), a list put in place of the circuit's, a copy; with the
    reference node removed too; and a name absent raising as before."""
    import copy

    c = _Two()
    c.append_node(Node('x'))
    _all_agree(c, extra=('nope',))
    c.connect_terminals(a=Node('n1'))
    ## (the replacing name first: a map kept past the edit lacks it, and a
    ## moved node looked up before it would remake the map)
    assert c.get_node_index('n1') == _scan_index(c, 'n1')
    _all_agree(c, extra=('a', 'nope'))

    class Pair(SubCircuit):
        terminals = ('p', 'm')

        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.add_node('mid')
            self['r1'] = R('p', 'mid', r=1.0)
            self['r2'] = R('mid', 'm', r=1.0)
    d = SubCircuit()
    d.add_nodes('in', 'out')
    d['x1'] = Pair('in', 'out')
    d['rg'] = R('out', gnd, r=1.0)
    _all_agree(d)
    assert d.get_node_index('out', gnd) == _scan_index(d, 'out') - (_scan_index(d, 'out') > _scan_index(d, gnd.name))
    del d['x1']
    _all_agree(d, extra=('x1.mid',))
    e = copy.copy(d)
    e.append_node(Node('late'))
    _all_agree(e)
    _all_agree(d, extra=('late',))
    d.nodes = [Node('q'), Node('r'), Node('s')]
    _all_agree(d, extra=('in',))
    p = SubCircuit()
    p.add_nodes('u', 'v', 'w')
    p.get_node_index('w')
    p.nodes.insert(0, p.nodes.pop())
    _all_agree(p)
