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
