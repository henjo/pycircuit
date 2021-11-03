# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test circuit module
"""

from nose.tools import *
import pycircuit.circuit.circuit 
from pycircuit.circuit.circuit import *
from pycircuit.circuit.elements import *
from pycircuit.circuit import AC, symbolic


from sympy import var, Symbol, simplify
import sympy
import numpy as np
from numpy.testing import * #assert_array_almost_equal, assert_array_equal
from copy import copy

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

    assert_equal(str(C(1, 0, gnd, c=sympy.Symbol('c'))),
                 "C('plus','minus',c=c)")

def test_print_netlist():
    """Test printing of netlist"""
    pycircuit.circuit.circuit.default_toolkit = numeric

    subc = generate_testcircuit()

    netlist = subc.netlist()
    print netlist
    
    refnetlist = \
""".subckt MySubC p m
  V1 internal gnd! VS v=0 vac=1 phase=0 noisePSD=0
  R1 p internal R r=1000.0 noisy=True
  R2 internal m R r=1000.0 noisy=True
  R3 internal gnd! R r=1000.0 noisy=True
.ends
I1 plus minus MySubC 
R1 plus minus R r=2000.0 noisy=True
R3 plus plus R r=1000.0 noisy=True"""

    assert_equal(netlist, refnetlist)    

def test_subcircuit_nodes():
    """Test node consistency of hierarchical circuit"""
    
    subc = generate_testcircuit()

    ## Check nodes of subc
    assert_equal(set(subc.nodes), 
                 set([Node('plus'), Node('minus'), Node('I1.internal'), 
                      gnd]))

    ## Check local names of subc
    assert_equal(subc.nodenames,
                 {'plus': Node('plus'), 'minus': Node('minus'),
                  'I1.internal': Node('I1.internal'),
                  'gnd': gnd})

    ## Check branches of subc
    assert_equal(subc.branches,
                 [Branch(Node('I1.internal'), gnd)])

    ## Check local names of I1
    assert_equal(subc['I1'].nodenames,
                 {'p': Node('p'), 'm': Node('m'),
                  'internal': Node('internal'),
                  'gnd': gnd})

    ## Check branches of I1
    assert_equal(subc['I1'].branches,
                 [Branch(Node('internal'), gnd)])

    ## Check nodes of I1
    assert_equal(set(subc['I1'].nodes), 
                 set([Node('p'), Node('m'), Node('internal'), 
                      gnd]))

    ## Check that first nodes of I1 are terminal nodes 
    assert_equal(subc['I1'].nodes[0:2], [Node('p'), Node('m')])

    ## Check terminal map
    assert_equal(subc.term_node_map['I1'], {'p':Node('plus'), 'm':Node('minus')})

    ## delete I1
    del subc['I1']
    
    ## Check nodes of subc
    assert_equal(set(subc.nodes), 
                 set([Node('plus'), Node('minus')]))

    ## Check local names of subc
    assert_equal(subc.nodenames,
                 {'plus': Node('plus'), 'minus': Node('minus')})
    
    ## Check terminal map
    assert_false('I1' in subc.term_node_map)

    ## Check nodes of R3
    assert_equal(subc['R3'].nodes,
                 [Node('plus'), Node('minus')])

def test_subcircuit_get_instance():
    cir = generate_testcircuit()

    assert_equal(cir[''], cir)
    assert_equal(cir['R1'], R('plus', 'minus', r=2e3))
    assert_equal(cir['I1.R1'], R('plus', 'minus', r=1e3))
    assert_raises(KeyError, lambda: cir['R10'])
    assert_raises(KeyError, lambda: cir['I1.R10'])
    assert_raises(KeyError, lambda: cir['I2.R10'])

def test_subcircuit_add_nodes_implicitly():
    subc = SubCircuit()

    ## Test to add nodes implicitly using node objects
    subc['R1'] = R(Node('a'), Node('b'))
    
    ## Check nodes of subc
    assert_equal(set(subc.nodes), 
                 set([Node('a'), Node('b')]))

    ## Check local names of subc
    assert_equal(subc.nodenames,
                 {'a': Node('a'), 'b': Node('b') })

    ## Test to add nodes implicitly using strings
    subc['R2'] = R('a', 'c')
    subc['R3'] = R('b', 1)
    
    ## Check nodes of subc
    assert_equal(set(subc.nodes), 
                 set([Node('a'), Node('b'), Node('c'), Node('1')]))

    ## Check local names of subc
    assert_equal(subc.nodenames,
                 {'a': Node('a'), 'b': Node('b'), 'c': Node('c'), 
                  '1': Node('1')})
    
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

    assert_equal(sympy.simplify(res.i('I1.plus')), (2 + C2*R3*s)/(1 + C2*R3*s))

    assert_equal(sympy.simplify(res.i('C2.plus')), s*R3*C2 / (1 + s*R3*C2))

            
def test_current_probing_wo_branch():
    """Test current probing with a current divider circuit without current probe"""

    s = sympy.Symbol('s')

    R1, C2, R3 = sympy.symbols(('R1', 'C2', 'R3'))

    cir = create_current_divider(R1,R3,C2)

    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)
    
    assert_equal(sympy.simplify(res.i('I1.plus')), (2 + C2*R3*s)/(1 + C2*R3*s))

    assert_equal(sympy.simplify(res.i('C2.plus')), s*R3*C2 / (1 + s*R3*C2))

def test_adddel_subcircuit_element():
    """add subcircuit element that contains a branch then delete it"""
    cir = SubCircuit()

    n1, = cir.add_nodes('n1')
    
    cir['R1'] = R(n1, gnd, r=1e3)
    
    cir['V'] = VS(n1, gnd)
    
    del cir['V']
    
    assert_equal(cir.elements.values(), [cir['R1']])
    assert_equal(cir.nodes, [n1,gnd])
    assert_equal(cir.branches, [])

def test_short_resistor():
    """Test shorting of instance terminals"""
    cir = SubCircuit()

    cir['R1'] = R(gnd, gnd)
    
    assert_equal(cir.G(np.zeros(1)), np.array([0]))
    
def test_copy_circuit():
    """Test to make a copy of circuit"""

    cir = generate_testcircuit()
    
    cir_copy = copy(cir)

    assert_equal(cir, cir_copy)

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

    print CircuitProxy(cir['I1'], cir, 'I1').terminalhook
    cir['I1'] = CircuitProxy(cir['I1'], cir, 'I1')
    
    assert_equal(cir['I1'].terminals, refcir['I1'].terminals)
    assert_equal(cir['I1'].non_terminal_nodes(), refcir['I1'].non_terminal_nodes())
    assert_equal(cir.nodes, refcir.nodes)
    assert_equal(cir.branches, refcir.branches)
    assert_equal(cir.n, refcir.n)

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

    assert_equal(a['R1'].iparv.r, 30)
    assert_equal(a['R1'].iparv.noisy, True)

    ## test 2 levels of hierarchy
    a['I1'] = A(x='x')
    a['I1']['R1'] = R(1,0, r='x+20')
    
    a.ipar.x = 30

    assert_equal(a['R1'].iparv.r, 40)
    assert_equal(a['I1']['R1'].iparv.r, 50)

def test_parameter_propagation_at_instantiation():
    """Test instance parameter value propagation through hierarchy at instantiation"""
    pycircuit.circuit.circuit.default_toolkit = symbolic

    class A(SubCircuit):
        instparams = [Parameter(name='resistance_value', desc = 'resistance value', default = 20)]

    a = A()
    
    a['R1'] = R(1,0, r='resistance_value + 10')

    assert_equal(a['R1'].iparv.r, 30)

    ## Verify that global parameters has lower priority than local parameters
    gp = ParameterDict(Parameter('resistance_value'))
    gp.resistance_value = 100
    a['R1'].update_iparv(a.iparv, gp)
    assert_equal(a['R1'].iparv.r, 30)

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

    assert_equal(b['R1'].iparv.r, 30)
    assert_equal(b['R2'].iparv.r, 40)
    
def test_global_parameters():
    a = SubCircuit()
    
    a['R1'] = R(1,0, r='R+10')
    
    globalparams = ParameterDict(Parameter('R'))

    globalparams.R = 20

    a.update_iparv(globalparams=globalparams)

    assert_equal(a['R1'].iparv.r, 30)

def test_replace_element():
    """Test node list consitency when replacing an element"""
    c = SubCircuit()
    c['VS'] = VS(1, gnd)
    assert_equal(set(c.nodes), set([Node('1'), gnd]))
    c['VS'] = VS(1, 0)
    assert_equal(set(c.nodes), set([Node('1'), Node('0')]))

def test_add_terminals():
    cir = SubCircuit()

    plus, common , minus = cir.add_nodes('plus', 'common', 'minus')
    
    cir.add_terminals(['plus','minus'])
    
    # The terminals list contains the names of all terminals
    assert_equal(cir.terminals,['plus','minus'])

    # The first k, k equal to the number of terminals, elements 
    # of the nodes list are the nodes connected to a terminal 
    node_list = []
    for node in cir.nodes:
        node_list.append(node.name)
    
    assert_equal(node_list,['plus','minus','common'])

def test_get_node():
    """Test  get_node""" 
    c = SubCircuit()
    out = c.add_node('out')
    c['V1'] = VS(out, gnd)
    assert_equal(c.get_node('V1.plus'), out)
