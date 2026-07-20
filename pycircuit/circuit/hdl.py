# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit import circuit
import pycircuit.utilities.param as param

import sympy
import sympy.printing.lambdarepr
import numpy as np

import inspect
from copy import copy

class Node(circuit.Node):
    @property
    def V(self):
        return Quantity('V', self)


class Branch(circuit.Branch):
    @property
    def V(self):
        return Quantity('V', self)

    @property
    def I(self):
        return Quantity('I', self)

class Parameter(param.Parameter, sympy.Symbol):
    pass

class ddt(sympy.Function):
    """Time derivative, d(x)/dt"""
    pass

class Quantity(sympy.AtomicExpr):
    """Reference to a voltage or current of a branch or node.

    Implemented as a sympy atom so it can take part in symbolic expressions
    (arithmetic, .subs(), Matrix.jacobian(), ...).
    """
    def __new__(cls, quantity, branch_or_node):
        if quantity not in ('V', 'I'):
            raise ValueError("quantity must be either 'V' or 'I'")
        if not isinstance(branch_or_node, (Node, Branch)):
            raise ValueError('branch_or_node must be a Branch or Node object')
        if quantity == 'I' and isinstance(branch_or_node, Node):
            raise ValueError('Current can only be taken on branches')

        obj = super().__new__(cls)
        obj.quantity = quantity
        obj.branch_or_node = branch_or_node
        return obj

    def _hashable_content(self):
        return (self.quantity, self.branch_or_node)

    @property
    def isnode(self): return isinstance(self.branch_or_node, Node)

    @property
    def isbranch(self): return isinstance(self.branch_or_node, Branch)

    def __repr__(self):
        if self.isbranch:
            return '%s(%s,%s)' % (self.quantity, self.branch_or_node.plus.name,
                                  self.branch_or_node.minus.name)
        return '%s(%s)' % (self.quantity, self.branch_or_node.name)

    __str__ = __repr__

class Statement():
    pass

class Contribution(Statement):
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = sympy.sympify(rhs)

    def nodes(self):
        """Return set of node objects referred to in lhs and rhs

        >>> a, b = Node('a'), Node('b')
        >>> b = Branch(a,b)
        >>> Contribution(b.I, 1e-3 * b.V).nodes()
        set([Node('a'), Node('b')])

        """
        
        nodes = set()

        for atom in self.lhs.atoms() | self.rhs.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch:
                    nodes.add(atom.branch_or_node.plus)
                    nodes.add(atom.branch_or_node.minus)
                else:
                    nodes.add(atom.branch_or_node)

        return nodes
        
    def contributions(self):
        """Return list of (node, iexpression, uexpression) tuples
        
        >>> a, b = Node('a'), Node('b')
        >>> b = Branch(a,b)
        >>> Contribution(b.I, 1e-3 * b.V).contributions()
        ((Node(a), 1e-3 * a.V - 1e-3 * b.V, 0), 
         (Node(b), -1e-3 * a.V + 1e-3 * b.V), 0)
        """
        if not isinstance(self.lhs, Quantity):
            raise ValueError('lhs must be a Quantity')

        rhs = self.rhs

        ## Split voltage of branches to voltages of nodes
        substdict = {}
        for atom in rhs.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch and atom.quantity == 'V':
                    branch = atom.branch_or_node
                    v = Quantity('V', branch.plus) - Quantity('V', branch.minus)
                    substdict[atom] = v
                    
        rhs = rhs.subs(substdict)
         
        ## Split i and u terms
        rhs = rhs.expand()

        if rhs.is_Add:
            terms = rhs.args
        else:
            terms = (rhs,)

        iterms = []
        uterms = []
        qterms = []
        for term in terms:
            if isconstant(term):
                uterms.append(term)
            elif isinstance(term, ddt):
                qterms.append(term.args[0])
            else:
                iterms.append(term)
                
        ## re-join terms 
        irhs = sympy.Add(*iterms)
        urhs = sympy.Add(*uterms)
        qrhs = sympy.Add(*qterms)
           
        if self.lhs.quantity == 'I':
            if self.lhs.isbranch:
                branch = self.lhs.branch_or_node
                return ((branch.plus, irhs, qrhs, urhs),
                        (branch.minus, -irhs, -qrhs, -urhs))

def _bind_x_method(func, paramnames):
    """Wrap a lambdified f(x, *params) as an instance method self.<name>(x)."""
    def method(self, x):
        return func(x, *[getattr(self.ipar, name) for name in paramnames])
    return method

def _bind_t_method(func, paramnames):
    """Wrap a lambdified f(t, *params) as an instance method self.u(t)."""
    def method(self, t=0.0):
        return func(t, *[getattr(self.ipar, name) for name in paramnames])
    return method

def generate_code(cls):
    """Return terminal names and lambdified i,u,q,G,C,CY functions.

    The analog() behaviour is turned into symbolic i/u/q vectors and the
    G/C Jacobians, which are then compiled to fast numpy functions with
    ``sympy.lambdify(..., cse=True)``.  The instance parameters are the
    trailing arguments of each compiled function and are supplied from
    ``self.ipar`` by the wrappers created in the metaclass.
    """

    ## Get arguments (terminals)
    terminalnames = inspect.getfullargspec(cls.analog)[0]

    ## Create node objects of the terminals
    terminalnodes = [Node(terminal) for terminal in terminalnames]

    ## Make a copy of analog method
    analogfunc = copy(cls.analog)

    ## Inject parameters into function globals so the analog body can refer to
    ## instance parameters by name.  They become plain symbols here and are
    ## bound to self.ipar.<name> values by the compiled-function wrappers.
    paramnames = [param.name for param in cls.instparams]
    paramsyms = [sympy.Symbol(name) for name in paramnames]
    analogfunc.__globals__.update(dict(zip(paramnames, paramsyms)))

    ## Call analog function
    statements = analogfunc(*terminalnodes)

    ## Create vector of current expressions for each node
    nodes = set()
    icontribs = {}
    ucontribs = {}
    qcontribs = {}
    for statement in statements:
        for node, icontrib, qcontrib, ucontrib in statement.contributions():
           if node in icontribs:
               icontribs[node] += icontrib
               ucontribs[node] += ucontrib
               qcontribs[node] += qcontrib
           else:
               icontribs[node] = icontrib
               ucontribs[node] = ucontrib
               qcontribs[node] = qcontrib

        nodes.update(statement.nodes())

    internalnodes = list(nodes - set(terminalnodes))

    nodes = terminalnodes + internalnodes

    ## Create a substitution dictionary that maps node voltages to symbols
    xvector = [sympy.Symbol('x_%d' % i) for i in range(len(nodes))]
    substdict = [(node.V, xsym) for node, xsym in zip(nodes, xvector)]

    ## Create i, u and q vectors
    ivector = [icontribs[node].subs(substdict) for node in nodes]
    qvector = [qcontribs[node].subs(substdict) for node in nodes]
    uvector = [ucontribs[node] for node in nodes]

    ## Calculate G as Jacobian of i
    icolvector = sympy.Matrix(ivector).T
    G = icolvector.jacobian(xvector)

    ## Calculate C as Jacobian matrix of q
    qcolvector = sympy.Matrix(qvector).T
    C = qcolvector.jacobian(xvector)

    CY = sympy.zeros(len(xvector))

    ## Compile to numpy functions.  A DeferredVector lets the state vector be
    ## passed as a single indexable argument, x, while the instance parameters
    ## follow as scalar arguments.  cse=True shares common subexpressions.
    x = sympy.DeferredVector('x')
    xsubs = dict(zip(xvector, [x[i] for i in range(len(xvector))]))

    def compile_x(expr):
        expr = expr.subs(xsubs) if hasattr(expr, 'subs') else \
            [e.subs(xsubs) for e in expr]
        return sympy.lambdify([x] + paramsyms, expr, modules='numpy', cse=True)

    t = sympy.Symbol('t')

    funcs = dict(
        i=compile_x(ivector),
        q=compile_x(qvector),
        G=compile_x(G),
        C=compile_x(C),
        CY=compile_x(CY),
        u=sympy.lambdify([t] + paramsyms, uvector, modules='numpy', cse=True),
    )

    return terminalnames, paramnames, funcs

class BehaviouralMeta(type):
    def __init__(cls, name, bases, dct):
        if 'analog' in dct:
            terminalnames, paramnames, funcs = generate_code(cls)

            ## Bind the compiled functions as instance methods, supplying the
            ## instance parameters from self.ipar.
            for methodname in ('i', 'q', 'G', 'C', 'CY'):
                setattr(cls, methodname,
                        _bind_x_method(funcs[methodname], paramnames))
            setattr(cls, 'u', _bind_t_method(funcs['u'], paramnames))

            ## Add terminals
            cls.terminals = terminalnames

                    
class Behavioural(circuit.Circuit, metaclass=BehaviouralMeta):
    """
    Behavioral circuit model

    The Behavioural is an extension of the Circuit class where an analogoue 
    circuit can be modelled at an abstract level that is similair to Verilog-A.
    
    The circuit behaviour is defined by the analog() method whose arguments
    are the terminal names and the voltages and currents are defined
    by calls to the contrib method. 
    
    Example
        class MyResistor(Behavioural):
            instparams = [param.Parameter(name='r', 
                          desc='Resistance', unit='ohm')]

            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, 1/self['r'] * b.I)

    The i(), q(), u(), G() and C() methods are then automatically generated
    from symbolic analysis of the expressions given to the contrib method.

    Using the hdl decorator on the definition of analog() adds some syntactic
    sugar that makes the behavior definition look more like Verilog-A.
    
    Example
        class MyResistor(Behavioural):
            instparams = [param.Parameter(name='r', 
                          desc='Resistance', unit='ohm')]

            @hdl
            def analog(plus, minus):
                I(plus, minus) <= 1/self['r'] * I(plus, minus)    
            
    """
    pass

def isconstant(expr):
    for atom in expr.atoms():
        if isinstance(atom, Quantity):
            return False
    return True
    
class Resistor(Behavioural):
      instparams = [Parameter(name='r', desc='Resistance', unit='ohm')]
      @staticmethod
      def analog(plus, minus):
          b = Branch(plus, minus)
          return Contribution(b.I, 1/r * b.V + 1),
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
