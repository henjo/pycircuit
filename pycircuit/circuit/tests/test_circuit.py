import pytest
# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test circuit module
"""

import pycircuit.circuit.circuit 
from pycircuit.circuit.circuit import *
from pycircuit.circuit.elements import *
from pycircuit.circuit import AC, symbolic


from sympy import var, Symbol, simplify
import sympy
import numpy as np
from numpy.testing import * #assert_array_almost_equal, assert_array_equal
from copy import copy

## (these tests set `circuit.default_toolkit` to the symbolic toolkit; the
## fixture puts it back -- tests/conftest.py)
pytestmark = pytest.mark.usefixtures('restore_default_toolkit')


def generate_testcircuit():
    subc = SubCircuit()
    plus, minus = subc.add_nodes('plus', 'minus')
    subc['R1'] = R(plus, minus, r=2e3)

    subc['R3'] = R(plus, plus)

    class MySubC(SubCircuit):
        terminals = ['p', 'm']

        def __init__(self, *args, **kvargs):
            super(MySubC, self).__init__(*args, **kvargs)
            
            internal = self.add_node('internal')
            
            self['R1'] = R(self.nodenames['p'], internal)
            self['R2'] = R(internal, self.nodenames['m'])
            self['V1'] = VS(internal, gnd)
            self['R3'] = R(internal, gnd)

    subc['I1'] = MySubC(plus, minus)

    return subc

def test_parallel():
    pycircuit.circuit.circuit.default_toolkit = numeric

    cir=SubCircuit()

    res = 1e3
    cir['R1'] = R(1, 2, res)
    cir['R2'] = R(1, 2, res)

    G = cir.G(np.array([0,0]))

    assert_array_equal(G, np.array([[2/res, -2/res],
                                    [-2/res, 2/res]]))

def test_print_element():
    pycircuit.circuit.circuit.default_toolkit = symbolic

    ## STAGE 10.3 added `ic` to C, so the repr carries it. Kept as an
    ## explicit expectation rather than made pattern-tolerant: the repr is
    ## meant to reconstruct the element, so a new parameter appearing in it is
    ## correct and worth pinning.
    assert str(C(1, 0, gnd, c=sympy.Symbol('c'))) == "C('plus','minus',c=c,ic=None)"

def test_print_netlist():
    """Test printing of netlist"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    subc = generate_testcircuit()

    netlist = subc.netlist()
    print(netlist)
    
    refnetlist = \
""".subckt MySubC p m
  R1 p internal R r=1000.0 noisy=True
  R2 internal m R r=1000.0 noisy=True
  V1 internal gnd! VS v=0 vac=1 phase=0 noisePSD=0
  R3 internal gnd! R r=1000.0 noisy=True
.ends
R1 plus minus R r=2000.0 noisy=True
R3 plus plus R r=1000.0 noisy=True
I1 plus minus MySubC """

    assert netlist == refnetlist    

def test_subcircuit_nodes():
    """Test node consistency of hierarchical circuit"""
    
    subc = generate_testcircuit()

    ## Check nodes of subc
    assert set(subc.nodes) == set([Node('plus'), Node('minus'), Node('I1.internal'), 
                      gnd])

    ## Check local names of subc
    assert subc.nodenames == {'plus': Node('plus'), 'minus': Node('minus'),
                  'I1.internal': Node('I1.internal'),
                  'gnd': gnd}

    ## Check branches of subc
    assert subc.branches == [Branch(Node('I1.internal'), gnd)]

    ## Check local names of I1
    assert subc['I1'].nodenames == {'p': Node('p'), 'm': Node('m'),
                  'internal': Node('internal'),
                  'gnd': gnd}

    ## Check branches of I1
    assert subc['I1'].branches == [Branch(Node('internal'), gnd)]

    ## Check nodes of I1
    assert set(subc['I1'].nodes) == set([Node('p'), Node('m'), Node('internal'), 
                      gnd])

    ## Check that first nodes of I1 are terminal nodes 
    assert subc['I1'].nodes[0:2] == [Node('p'), Node('m')]

    ## Check terminal map
    assert subc.term_node_map['I1'] == {'p':Node('plus'), 'm':Node('minus')}

    ## delete I1
    del subc['I1']
    
    ## Check nodes of subc
    assert set(subc.nodes) == set([Node('plus'), Node('minus')])

    ## Check local names of subc
    assert subc.nodenames == {'plus': Node('plus'), 'minus': Node('minus')}
    
    ## Check terminal map
    assert not 'I1' in subc.term_node_map

    ## Check nodes of R3
    assert subc['R3'].nodes == [Node('plus'), Node('minus')]

def test_subcircuit_get_instance():
    cir = generate_testcircuit()

    assert cir[''] == cir
    assert cir['R1'] == R('plus', 'minus', r=2e3)
    assert cir['I1.R1'] == R('plus', 'minus', r=1e3)
    with pytest.raises(KeyError):
        cir['R10']
    with pytest.raises(KeyError):
        cir['I1.R10']
    with pytest.raises(KeyError):
        cir['I2.R10']

def test_subcircuit_add_nodes_implicitly():
    subc = SubCircuit()

    ## Test to add nodes implicitly using node objects
    subc['R1'] = R(Node('a'), Node('b'))
    
    ## Check nodes of subc
    assert set(subc.nodes) == set([Node('a'), Node('b')])

    ## Check local names of subc
    assert subc.nodenames == {'a': Node('a'), 'b': Node('b') }

    ## Test to add nodes implicitly using strings
    subc['R2'] = R('a', 'c')
    subc['R3'] = R('b', 1)
    
    ## Check nodes of subc
    assert set(subc.nodes) == set([Node('a'), Node('b'), Node('c'), Node('1')])

    ## Check local names of subc
    assert subc.nodenames == {'a': Node('a'), 'b': Node('b'), 'c': Node('c'), 
                  '1': Node('1')}
    
def create_current_divider(R1,R3,C2):
    cir = SubCircuit()

    n1,n2 = cir.add_nodes('n1', 'n2')
    
    class MySubC(SubCircuit):
        terminals = ['plus', 'minus']

        def __init__(self, *args, **kvargs):
            super(MySubC, self).__init__(*args, **kvargs)

            self['R3'] = R(self.nodenames['plus'], self.nodenames['minus'], r=R3)
            self['I2'] = IS(self.nodenames['plus'], self.nodenames['minus'], iac=1)


    cir['IS'] = IS(gnd,n1, iac=2)
    cir['R1'] = R(n1, n2, r=R1)
    cir['I1'] = MySubC(n2, gnd)
    cir['C2'] = C(n2, gnd, c=C2)
 
    return cir

def test_current_probing():
    """Test current probing with a current divider circuit"""
    pycircuit.circuit.circuit.default_toolkit = symbolic
    
    s = sympy.Symbol('s')

    R1, R3, C2 = sympy.symbols(('R1', 'R3', 'C2'))

    cir = create_current_divider(R1,R3,C2)
    
    cir = cir.save_current('I1.plus')
    
    assert cir.get_terminal_branch('I1.plus') is not None
    
    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert sympy.simplify(res.i('I1.plus')) == (2 + C2*R3*s)/(1 + C2*R3*s)

    assert sympy.simplify(res.i('C2.plus')) == s*R3*C2 / (1 + s*R3*C2)

            
def test_current_probing_wo_branch():
    """Test current probing with a current divider circuit without current probe"""
    import pycircuit.circuit.circuit
    pycircuit.circuit.circuit.default_toolkit = symbolic

    s = sympy.Symbol('s')

    R1, C2, R3 = sympy.symbols(('R1', 'C2', 'R3'))

    cir = create_current_divider(R1,R3,C2)

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)
    
    assert sympy.simplify(res.i('I1.plus')) == (2 + C2*R3*s)/(1 + C2*R3*s)

    assert sympy.simplify(res.i('C2.plus')) == s*R3*C2 / (1 + s*R3*C2)

def test_adddel_subcircuit_element():
    """add subcircuit element that contains a branch then delete it"""
    cir = SubCircuit()

    n1, = cir.add_nodes('n1')
    
    cir['R1'] = R(n1, gnd, r=1e3)
    
    cir['V'] = VS(n1, gnd)
    
    del cir['V']
    
    assert list(cir.elements.values()) == [cir['R1']]
    assert cir.nodes == [n1,gnd]
    assert cir.branches == []

def test_short_resistor():
    """Test shorting of instance terminals"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    cir = SubCircuit()

    cir['R1'] = R(gnd, gnd)

    ## G is the n-by-n MNA conductance matrix; for the single gnd node
    ## that is a 1-by-1 zero matrix.
    assert cir.G(np.zeros(1)) == np.zeros((1, 1))
    
def test_copy_circuit():
    """Test to make a copy of circuit"""

    cir = generate_testcircuit()
    
    cir_copy = copy(cir)

    assert cir == cir_copy

def test_VCCS_tied():
    """Test VCCS with some nodes tied together"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    cir = SubCircuit()

    n3,n2 = cir.add_nodes('3','2')

    gm1 = sympy.Symbol('gm1')

    cir['gm'] = VCCS(gnd, n3, n2, n3, gm = gm1)   
    
    assert_array_equal(cir.G(np.zeros(cir.n)),
                       np.array([[gm1, 0, -gm1],
                                 [-gm1, 0, gm1],
                                 [0, 0, 0]]))

    
def test_proxy():
    pycircuit.circuit.circuit.default_toolkit = symbolic
    
    refcir = generate_testcircuit()
    
    cir = generate_testcircuit()

    print(CircuitProxy(cir['I1'], cir, 'I1').terminalhook)
    cir['I1'] = CircuitProxy(cir['I1'], cir, 'I1')
    
    assert cir['I1'].terminals == refcir['I1'].terminals
    assert cir['I1'].non_terminal_nodes() == refcir['I1'].non_terminal_nodes()
    assert cir.nodes == refcir.nodes
    assert cir.branches == refcir.branches
    assert cir.n == refcir.n

    for method in ['G', 'C', 'i', 'q']:
        assert_array_equal(getattr(cir, method)(np.zeros(cir.n)),
                           getattr(refcir, method)(np.zeros(cir.n)),
                           )

    assert_array_equal(cir.CY(np.zeros(cir.n),1), refcir.CY(np.zeros(cir.n),1))

def test_parameter_propagation():
    """Test instance parameter value propagation through hierarchy"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    class A(SubCircuit):
        instparams = [Parameter('x')]

    a = A()

    a['R1'] = R(1,0, r='10+x')

    a.ipar.x = 20

    assert a['R1'].iparv.r == 30
    assert a['R1'].iparv.noisy == True

    ## test 2 levels of hierarchy
    a['I1'] = A(x='x')
    a['I1']['R1'] = R(1,0, r='x+20')
    
    a.ipar.x = 30

    assert a['R1'].iparv.r == 40
    assert a['I1']['R1'].iparv.r == 50

def test_parameter_propagation_at_instantiation():
    """Test instance parameter value propagation through hierarchy at instantiation"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    class A(SubCircuit):
        instparams = [Parameter(name='resistance_value', desc = 'resistance value', default = 20)]

    a = A()
    
    a['R1'] = R(1,0, r='resistance_value + 10')

    assert a['R1'].iparv.r == 30

    ## Verify that global parameters has lower priority than local parameters
    gp = ParameterDict(Parameter('resistance_value'))
    gp.resistance_value = 100
    a['R1'].update_iparv(a.iparv, gp)
    assert a['R1'].iparv.r == 30

def _resolved_hierarchy(nested=True):
    """A parameterised circuit `a` (x = 20) holding R1 (r = '10+x'), and
    (`nested`) an instance I1 of the same class (x = 'x') holding R1
    (r = 'x+20')."""
    pycircuit.circuit.circuit.default_toolkit = numeric

    class A(SubCircuit):
        instparams = [Parameter('x', default=0.0)]  # noqa: RUF012

    a = A(x=20.0)
    a['R1'] = R(1, 0, r='10+x')
    if nested:
        a['I1'] = A(x='x')
        a['I1']['R1'] = R(1, 0, r='x+20')
    return a


def test_an_element_s_own_write_resolves_against_its_parent():
    """An element's own parameter write re-resolves its parameters against
    the parent it was last resolved in: another parameter given as an
    expression of the parent's keeps its value, and a new expression
    written onto it resolves (until 2026-10-22 it resolved against nothing:
    r fell back to its default, 1 kohm, after `noisy` was written)."""
    a = _resolved_hierarchy()
    assert a['R1'].iparv.r == 30.0
    a['R1'].ipar.noisy = False
    assert a['R1'].iparv.r == 30.0
    a['R1'].ipar.r = '5+x'
    assert a['R1'].iparv.r == 25.0
    ## (a nested instance's own write: its parameter from the top's, its
    ## element's from it)
    a['I1'].ipar.x = 'x*2'
    assert a['I1'].iparv.x == 40.0
    assert a['I1']['R1'].iparv.r == 60.0
    a['I1']['R1'].ipar.noisy = False
    assert a['I1']['R1'].iparv.r == 60.0


def test_a_copied_circuit_s_elements_resolve_against_the_copy():
    """(A flat circuit: copying one that holds a subcircuit instance raises
    KeyError in `update_node_map`, before and after.)"""
    a = _resolved_hierarchy(nested=False)
    b = copy(a)
    b.ipar.x = 50.0
    assert b['R1'].iparv.r == 60.0
    b['R1'].ipar.noisy = False
    assert b['R1'].iparv.r == 60.0
    assert a['R1'].iparv.r == 30.0


def test_an_unchanged_resolution_context_is_not_rewritten():
    """The remembered context is stored only when it changes: a re-resolution
    in the same parent leaves the element's dict as it was (the dict-watch
    stamps read it)."""
    a = _resolved_hierarchy()
    ctx = a['R1'].__dict__['_resolved_in']
    a.update_iparv()
    a['R1'].ipar.noisy = False
    assert a['R1'].__dict__['_resolved_in'] is ctx


def _divider(added_first):
    """A source across a two-resistor divider in a subcircuit instance,
    the instance given its elements before or AFTER it is added."""
    from pycircuit.circuit.dcanalysis import DC  # noqa: F401
    pycircuit.circuit.circuit.default_toolkit = numeric

    class Div(SubCircuit):
        terminals = ('a', 'b')

    top = SubCircuit()
    n1 = top.add_node('in')
    top['V'] = VS(n1, gnd, v=2.0)
    d = Div(n1, gnd)
    if added_first:
        top['X'] = d
    d['R1'] = R('a', 'm', r=1e3)
    d['R2'] = R('m', 'b', r=3e3)
    if not added_first:
        top['X'] = d
    return top


def test_an_instance_given_its_elements_after_it_was_added():
    """The parent follows the instance's later elements -- its internal
    node and the answer are the elements-first circuit's (until 2026-10-22
    the parent kept its structure from the addition: DC and copy raised
    KeyError on Node('X.m'))."""
    from pycircuit.circuit.dcanalysis import DC
    late, early = _divider(True), _divider(False)
    assert late.nodes == early.nodes and late.branches == early.branches
    got = DC(late).solve()
    assert np.asarray(got.x).tobytes() == np.asarray(DC(early).solve().x).tobytes()
    assert np.isclose(float(got.v('X.m', gnd)), 1.5, rtol=1e-12)
    c = copy(late)
    assert np.asarray(DC(c).solve().x).tobytes() == np.asarray(got.x).tobytes()


def test_a_nested_instance_s_later_structure_reaches_the_top():
    """Two levels: an element added to the inner instance after both were
    added reaches the top circuit; one removed from it leaves it."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    top = SubCircuit()
    mid = SubCircuit()
    inner = SubCircuit()
    top['M'] = mid
    mid['I'] = inner
    inner['R'] = R('p', 'q', r=1.0)
    assert {'M.I.p', 'M.I.q'} <= {n.name for n in top.nodes}
    del inner['R']
    assert not {'M.I.p', 'M.I.q'} & {n.name for n in top.nodes}


def test_a_deleted_instance_leaves_the_later_branch_rows_right():
    """After a deletion each later instance's branch rows are its own (the
    spans were left unshifted until 2026-10-22: a row past the end), and
    of two parallel inductors the deleted one's branch goes."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    a = c.add_node('a')
    c['V0'] = VS(a, gnd, v=1.0)
    c['L1'] = L(a, gnd, L=1e-6)
    c['L2'] = L(a, gnd, L=2e-6)
    del c['V0']
    nn = len(c.nodes)
    assert c.instance_branch_indices('L1') == [nn]
    assert c.instance_branch_indices('L2') == [nn + 1]
    assert 'V0' not in c._instance_branch_span
    del c['L1']
    assert c.instance_branch_indices('L2') == [nn]
    assert len(c.branches) == 1


def test_parameter_at_instantiation_with_add_instance ():
    """Test instance parameter value propagation through hierarchy with add_instance"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    a = SubCircuit()

    a['R1'] = R(1,0, r=10)

    # copy 'R1' resistor
    res  = copy(a['R1'])
    res2 = copy(a['R1'])

    # change resistance value
    res.ipar.set(r=30)
    res2.ipar.set(r= 40)

    b = SubCircuit()
    # insert copied instance into circuit b
    b.add_instance('R1',res, plus='plus', minus='minus')
    b.add_instance('R2',res2, plus='plus', minus='minus')

    assert b['R1'].iparv.r == 30
    assert b['R2'].iparv.r == 40
    
def test_global_parameters():
    a = SubCircuit()
    
    a['R1'] = R(1,0, r='R+10')
    
    globalparams = ParameterDict(Parameter('R'))

    globalparams.R = 20

    a.update_iparv(globalparams=globalparams)

    assert a['R1'].iparv.r == 30

def test_replace_element():
    """Test node list consitency when replacing an element"""
    c = SubCircuit()
    c['VS'] = VS(1, gnd)
    assert set(c.nodes) == set([Node('1'), gnd])
    c['VS'] = VS(1, 0)
    assert set(c.nodes) == set([Node('1'), Node('0')])

def test_add_terminals():
    cir = SubCircuit()

    plus, common , minus = cir.add_nodes('plus', 'common', 'minus')
    
    cir.add_terminals(['plus','minus'])
    
    # The terminals list contains the names of all terminals
    assert cir.terminals == ['plus','minus']

    # The first k, k equal to the number of terminals, elements 
    # of the nodes list are the nodes connected to a terminal 
    node_list = []
    for node in cir.nodes:
        node_list.append(node.name)
    
    assert node_list == ['plus','minus','common']

def test_get_node():
    """Test  get_node""" 
    c = SubCircuit()
    out = c.add_node('out')
    c['V1'] = VS(out, gnd)
    assert c.get_node('V1.plus') == out


def test_extract_i_at_an_independent_source_terminal_is_the_source_current():
    """2026-09-08.  `extract_i` derives the terminal current as
    `i(x) + u(t) + C xdot` at the terminal node -- applied at ELEMENT level,
    since `SubCircuit.extract_i` hands the terminal to the element with its
    own slice of x.  The time-domain branches used to call `u(t)` bare,
    which returns ZEROS, so an independent source's own terminal read 0.0
    while every other element was unaffected.  Pinned both ways: the
    resistor's terminal (unchanged) and the source's (was 0.0).
    """
    import numpy as np
    from pycircuit.circuit import circuit
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import IS, R
    from pycircuit.circuit.dcanalysis import DC
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('a')
    cir['I1'] = IS('a', gnd, i=2e-3)
    cir['R1'] = R('a', gnd, r=1e3)
    x = np.asarray(DC(cir).solve().x, dtype=float)
    assert abs(float(cir.extract_i(x, 'R1.plus')) - (-2e-3)) < 1e-12
    assert abs(float(cir.extract_i(x, 'I1.plus')) - (+2e-3)) < 1e-12, \
        'the source terminal read 0.0 before the analysis flag was passed'
    assert abs(float(cir.extract_i(x, 'I1.minus')) - (-2e-3)) < 1e-12


def test_circuit_proxy_forwards_every_argument():
    """`CircuitProxy.u` raised NameError on every call (`self.device.u(x, epar)`
    with no `x` in scope) -- latent, since no shipped path calls a proxy's
    `u` -- and the obvious `x -> t` fix would have made the wrapped SOURCE
    VANISH silently by dropping `analysis` (docs session, 2026-09-09).  `CY`
    passed `epar` where the device expects `w`.  So: the proxy's `u` with
    `analysis='dc'` equals the device's, its bare `u` equals the device's
    bare `u` (zeros, by the classical gate), and `CY` forwards `w`.
    """
    from pycircuit.circuit import circuit
    from pycircuit.circuit.circuit import CircuitProxy, SubCircuit, gnd
    from pycircuit.circuit.elements import R, IS
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('a')
    cir['I1'] = IS('a', gnd, i=2e-3, noisePSD=1e-6)
    cir['R1'] = R('a', gnd, r=1e3)
    dev = cir['I1']
    proxy = CircuitProxy(dev, cir, 'I1')
    x = np.zeros(dev.n)
    assert np.allclose(np.asarray(proxy.u(0.0, analysis='dc'), float),
                       np.asarray(dev.u(0.0, analysis='dc'), float))
    assert np.any(np.asarray(dev.u(0.0, analysis='dc'), float) != 0.0)
    assert np.allclose(np.asarray(proxy.u(0.0), float), np.asarray(dev.u(0.0), float))
    for w in (1.0, 1e6):
        assert np.allclose(np.asarray(proxy.CY(x, w), complex), np.asarray(dev.CY(x, w), complex))
    assert np.allclose(np.asarray(proxy.G(x), float), np.asarray(dev.G(x), float))
    assert np.allclose(np.asarray(proxy.i(x), float), np.asarray(dev.i(x), float))

