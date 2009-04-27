# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

""" Test circuit module
"""

from nose.tools import *
from pycircuit.circuit.circuit import *
from pycircuit.circuit.elements import *
from pycircuit.circuit import AC, symbolic

import sympy
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from numpy.testing.decorators import slow
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
    cir=SubCircuit()

    res = 1e3
    cir['R1'] = R(1, 2, res)
    cir['R2'] = R(1, 2, res)

    G = cir.G(array([0,0]))

    assert_array_equal(G, array([[2/res, -2/res],
                                 [-2/res, 2/res]]))

def test_print_element():
    assert_equal(str(C(1, 0, gnd, c=sympy.Symbol('c'))),
                 "C('plus','minus',c=c)")

def test_print_netlist():
    """Test printing of netlist"""

    subc = generate_testcircuit()

    netlist = subc.netlist()
    print netlist
    
    refnetlist = \
""".subckt MySubC p m
  V1 internal gnd! VS v=1 vac=1 noisePSD=0
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
    
    s = sympy.Symbol('s')

    R1, R3, C2 = sympy.symbols('R1', 'R3', 'C2')

    cir = create_current_divider(R1,R3,C2)
    
    cir = cir.save_current('I1.plus')
    
    assert cir.get_terminal_branch('I1.plus') != None
    
    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)

    assert_equal(sympy.simplify(res.i('I1.plus')), (2 + C2*R3*s)/(1 + C2*R3*s))

    assert_equal(sympy.simplify(res.i('C2.plus')), s*R3*C2 / (1 + s*R3*C2))

            
def test_current_probing_wo_branch():
    """Test current probing with a current divider circuit without current probe"""

    s = sympy.Symbol('s')

    R1, C2, R3 = sympy.symbols('R1', 'C2', 'R3')

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
    
    assert_equal(cir.G(zeros(1)), array([0]))
    
def test_copy_circuit():
    """Test to make a copy of circuit"""

    cir = generate_testcircuit()
    
    cir_copy = copy(cir)

    assert_equal(cir, cir_copy)

def test_VCCS_tied():
    """Test VCCS with some nodes tied together"""
    cir = SubCircuit()

    n3,n2 = cir.add_nodes('3','2')

    gm1 = sympy.Symbol('gm1')

    cir['gm'] = VCCS(gnd, n3, n2, n3, gm = gm1)   
    
    assert_array_equal(cir.G(zeros(cir.n)),
                       array([[gm1, 0, -gm1],
                              [-gm1, 0, gm1],
                              [0, 0, 0]]))

    
def test_proxy():
    
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
        assert_array_equal(getattr(cir, method)(zeros(cir.n)),
                           getattr(refcir, method)(zeros(cir.n)),
                           )

    assert_array_equal(cir.CY(zeros(cir.n),1), refcir.CY(zeros(cir.n),1))

def test_parameter_propagation():
    """Test instance parameter value propagation through hierarchy"""

    class A(SubCircuit):
        instparams = [Parameter('x')]

    a = A()

    a['R1'] = R(1,0, r=Parameter('x') + 10)

    a.ipar.x = 20

    assert_equal(a['R1'].ipar.r, 30)
    assert_equal(a['R1'].ipar.noisy, True)

    ## test 2 levels of hierarchy
    a['I1'] = A(x=Parameter('x'))
    a['I1']['R1'] = R(1,0, r=Parameter('x') + 20)
    
    a.ipar.x = 30

    assert_equal(a['R1'].ipar.r, 40)
    assert_equal(a['I1']['R1'].ipar.r, 50)
    
def test_design_variables():
    a = SubCircuit()
    
    a['R1'] = R(1,0, r=Variable('R')+10)
    
    ipars = ParameterDict()
    variables = ParameterDict(Variable('R'))

    variables.R = 20

    a.update_ipar(ipars, variables)

    assert_equal(a['R1'].ipar.r, 30)

def test_replace_element():
    """Test node list consitency when replacing an element"""
    c = SubCircuit()
    c['VS'] = VS(1, gnd)
    assert_equal(set(c.nodes), set([Node('1'), gnd]))
    c['VS'] = VS(1, 0)
    assert_equal(set(c.nodes), set([Node('1'), Node('0')]))
    
