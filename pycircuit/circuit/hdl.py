import circuit
import pycircuit.utilities.param as param

import sympy
import sympy.printing.lambdarepr
import numpy as npy

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

class dtt(sympy.Function):
    """Time derivative, d(x)/dt"""
    pass

class Quantity(circuit.Quantity, sympy.basic.Atom):
    pass

class Statement(object):
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
            elif isinstance(term, dtt):
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

class NumpyPrinter(sympy.printing.StrPrinter):
    def _print_Matrix(self, expr):
        return "npy.array([%s])"%expr._format_str(self._print, ",")

def methodstr(name, args, expr):
    """Returns a string that can be evaluated to an instance method
    
    >>> x = sympy.Symbol('x')
    >>> methodstr('G', x, x**2 * sympy.Symbol('self.p'))
    'def G(self, x): return self.p*x**2'
    
    """

   # Transform everything to strings.
    expr = NumpyPrinter().doprint(expr)
    if isinstance(args, str):
        pass
    elif hasattr(args, "__iter__"):
        args = ",".join(str(a) for a in args)
    else:
        args = str(args)

    return "def %s(self, %s): return %s" % (name, args, expr)

def generate_code(cls):
    """Returns terminal names and i,u,q,G,C,CY method strings from class obj"""

    ## Get arguments (terminals)
    terminalnames = inspect.getargspec(cls.analog)[0]

    ## Create node objects of the terminals
    terminalnodes = [Node(terminal) for terminal in terminalnames]

    ## Make a copy of analog method
    analogfunc = copy(cls.analog)

    ## Inject parameters into function globals
    params = dict((param.name, sympy.Symbol('self.ipar.' + param.name))
                  for param in cls.instparams)
    analogfunc.func_globals.update(params)

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

    ## Create a substitution dictionary that maps node voltages 
    ## to symbols
    xvector = [sympy.Symbol('x[%d]'%i) for i in range(len(nodes))]
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

    ## Create Circuit methods
    ifuncstr = methodstr('i', 'x', ivector)
    ufuncstr = methodstr('u', 't', uvector)
    qfuncstr = methodstr('q', 'x', qvector)
    Gfuncstr = methodstr('G', 'x', G)
    Cfuncstr = methodstr('C', 'x', C)
    CYfuncstr = methodstr('CY', 'x', CY)

    return terminalnames, ifuncstr, ufuncstr, qfuncstr, Gfuncstr, Cfuncstr, \
        CYfuncstr
        
class BehaviouralMeta(type):
    def __init__(cls, name, bases, dct):
        if 'analog' in dct:
            ## Generate code for the Circuit methods as strings
            strings = generate_code(cls)
            terminalnames, ifuncstr, ufuncstr, qfuncstr, Gfuncstr, \
                Cfuncstr, CYfuncstr = strings

            ## Create methods
            methodnames = ('i', 'u', 'q', 'G', 'C', 'CY')
            for methodname, codestring in zip(methodnames, strings[1:]):
                funcdef_code = compile(codestring, '<stdin>', 'exec')
                namespace = {'npy': npy}
                eval(funcdef_code, namespace)
                setattr(cls, methodname, namespace[methodname])

            ## Add terminals
            cls.terminals = terminalnames

                    
class Behavioural(circuit.Circuit):
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
    
    __metaclass__ = BehaviouralMeta

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
