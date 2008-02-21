import unittest
from numpy import *
import pylab
import sympy
from sympy import Symbol, Matrix, symbols, simplify

class Node:
    """A node is a region in an electric circuit where the voltage is the same.
    
    A node can have a name
    """
    def __init__(self, name=None):
        self.name = name

    def __repr__(self):
        return self.__class__.__name__ + '(\'' + self.name + '\')'

class Branch:
    """A branch connects two nodes.
    
    A branch is used in modified nodal analysis to describe components that defines a voltage between
    two nodes as a function of current flowing between these nodes. Examples are voltage sources and inductors.
    Positive current through a branch is defined as a current flowing from plus to minus.a
    
    """
    def __init__(self, plus, minus, name=None):
        """Initiate a branch

        Arguments:
        plus -- Node object connected to the postive terminal of the branch
        minus -- Node object connected to the negative terminal of the branch

        Keyword arguments:
        name -- branch name

        """

        self.plus = plus
        self.minus = minus
        self.name = name

class Terminal:
    def __init__(self, name):
        """A Terminal connects the nodes in a Circuit to nodes in superior hierarchy levels"""
        self.name = name

class Circuit:
    """The circuit class describes a full circuit, subcircuit or a single component. 

    It contains a list of nodes,terminals and parameters.
    The terminals connect nodes inside the circuit to the nodes outside the circuit. When the circuit is
    instanciated, the outside nodes are passed to the object via the terminals.
       
    Attributes:
    nodes      -- A list that contains Node objects. If the circuit does not contain any internal nodes
                  the length is the same as the number of terminals.
    branches   -- A list of Branch objects. The solver will solve for the currents through the branches.
                  terminals  A list that contains Terminal objects
    parameters -- A dictionary of parameter values keyed by the parameter name
    nodenames  -- A dictionary that maps a local node name to the node object in nodes. If the node is
                  connnected to superior hierarchy levels through a terminal the terminal name must
                  be the same as the local node name


    >>> circuit=Circuit(terminals=["plus", "minus"])
    >>> circuit.terminals
    ['plus', 'minus']
    
    
    """
    def __init__(self, terminals=[]):
        self.nodes=[]
        self.terminals=[]
        self.nodenames={}
        self.branches=[]
        self.parameters={}

        map(self.addTerminal, terminals)
        
    def addNode(self, name=None):
        """Create an internal node in the circuit and return the new node

        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.nodes
        [Node('n1')]
        >>> 'n1' in c.nodenames
        True
        
        """
        newnode = Node(name)
        self.nodes.append(newnode)
        if name != None:
            self.nodenames[name] = newnode
        return newnode

    def getNode(self, name):
        """Find a node by name.
        
        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.getNode('n1')
        Node('n1')
        
        """
        return self.nodenames[name]

    def addTerminal(self, name):
        """Add a new terminal to the circuit and create a node connected to it

        >>> circuit=Circuit()
        >>> circuit.terminals
        []
        >>> circuit.addTerminal("plus")
        >>> circuit.terminals
        ['plus']
        >>> circuit.nodenames["plus"]
        Node('plus')
        
        """
        self.terminals.append(name)
        self.addNode(name)

    def connectNode(self, terminal, node):
        """Connect an external node to a terminal

        >>> circuit=Circuit(terminals=["plus", "minus"])
        >>> node=Node("outsidenode")
        >>> circuit.nodenames["plus"]
        Node('plus')
        >>> circuit.connectNode("plus", node)
        >>> circuit.nodenames["plus"]
        Node('outsidenode')
        
        """
        if node != None:
            if not terminal in self.terminals:
                raise ValueError('terminal '+str(terminal)+' is not defined')
            self.nodes.remove(self.nodenames[terminal])
            self.nodes.append(node)
            self.nodenames[terminal] = node

    def G(self, x):
        """Calculate the G ((trans)conductance) matrix of the circuit given the x-vector"""
        N=len(self.nodes)+len(self.branches)
        return zeros((N,N), dtype=object)

    def C(self, x):
        """Calculate the C ((trans)capacitance) matrix of the circuit given the x-vector"""
        N=len(self.nodes)+len(self.branches)
        return zeros((N,N), dtype=object)

    def U(self):
        """Calculate the U (circuit input) column-vector the circuit"""
        N=len(self.nodes)+len(self.branches)
        return zeros((N,1), dtype=object)

    def solvedc(self):
        N=len(self.nodes)+len(self.branches)

        G=self.G(zeros((N,1)))
        U=self.U()

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        ignd = self.nodes.index(gnd)
        G=removeRowCol(G, ignd)
        U=removeRowCol(U, ignd)

        return linalg.solve(G,-U)

    def solveac(self, freqs):
        N=len(self.nodes)+len(self.branches)

        G=self.G(zeros((N,1)))
        C=self.C(zeros((N,1)))
        U=self.U()

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        ignd = self.nodes.index(gnd)
        G=removeRowCol(G, ignd).astype(complex)
        C=removeRowCol(C, ignd).astype(complex)
        U=removeRowCol(U, ignd).astype(complex)

        out = []
        for f in freqs:
            out.append(linalg.solve(2j*pi*f*C + G, U))
        return out

    def solvesymbolic(self):
        N=len(self.nodes)+len(self.branches)

        G=self.G(zeros((N,1)))
        C=self.C(zeros((N,1)))
        U=self.U()

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        ignd = self.nodes.index(gnd)
        G=removeRowCol(G, ignd)
        C=removeRowCol(C, ignd)
        U=removeRowCol(U, ignd)

        G=sympy.Matrix(G)
        C=sympy.Matrix(C)
        U=sympy.Matrix(U)
        outputvariables = map(Symbol, ['v'+node.name for node in self.nodes if not node is gnd]+
                              ["I%d"%i for i in range(len(self.branches))])
        return sympy.solve_linear_system((Symbol('s')*C+G).row_join(-U), outputvariables)

class SubCircuit(Circuit):
    """
    SubCircuit is container for circuits.
    Attributes:
      elements          dictionary of Circuit objects keyed by its instance name
      elementnodemap    list of translations between node indices of the elements to the
                        node index in the SubCircuit object.
      elementbranchmap  list of translations between branch indices of the elements to the
                        branch index in the SubCircuit object.
    """
    def __init__(self):
        Circuit.__init__(self)
        self.elements = {}
        self.elementnodemap = []
        self.elementbranchmap = []
        
    def __setitem__(self, instancename, element):
        """Adds an instance to the circuit"""
        self.elements[instancename] = element
        
        # Add nodes and branches from new element
        for node in element.nodes:
            if not node in self.nodes:
                self.nodes.append(node)
        self.branches.extend(element.branches)

        self.updateNodeMap()

    def getNode(self, name):
        """Find a node by name.

        The name can be a hierachical name and the notation is '/I1/I2/net1' for
        a net net1 in instance I2 of instance I1
        
        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.getNode('n1')
        Node('n1')

        >>> c1 = SubCircuit()
        >>> c2 = SubCircuit()
        >>> c1['I1'] = c2
        >>> n1 = c2.addNode("net1")
        >>> c1.getNode('/I1/net1')
        Node('net1')
        
        """
        hierlevels = [part for part in name.split('/') if part != '']
            
        if len(hierlevels)==1:
            return self.nodenames[hierlevels[0]]
        else:
            return self.elements[hierlevels[0]].getNode('/'.join(['']+hierlevels[1:]))

    def updateNodeMap(self):
        """Update the elementnodemap attribute"""

        self.elementnodemap = {}
        for instance, element in self.elements.items():
            self.elementnodemap[instance] = [self.nodes.index(node) for node in element.nodes] + \
                                            [self.branches.index(branch)+len(self.nodes) for branch in element.branches]

    def G(self, x):
        N=len(self.nodes)+len(self.branches)
        G=zeros((N,N), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            G[[[i] for i in nodemap], nodemap] += element.G(x)
            
        return G

    def C(self, x):
        N=len(self.nodes)+len(self.branches)
        C=zeros((N,N), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            C[[[i] for i in nodemap], nodemap] += element.C(x)
        return C

    def U(self):
        N=len(self.nodes)+len(self.branches)
        U=zeros((N,1), dtype=object)

        self.updateNodeMap()
        
        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            U[nodemap] += element.U()
        return U

class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c.G(0)
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)
    
    """
    def __init__(self, plus, minus, r=1e3):
        Circuit.__init__(self, terminals=['plus', 'minus'])
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['r']=r
        
    def G(self, x):
        return array([[1/self.parameters['r'] , -1/self.parameters['r']],
                       [-1/self.parameters['r'], 1/self.parameters['r']] ])

class C(Circuit):
    """Capacitor"""
    def __init__(self, plus, minus, C=0.0):
        Circuit.__init__(self, terminals=['plus', 'minus'])
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['C']=C
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def C(self, x):
        return array([[self.parameters['C'] , -self.parameters['C']],
                       [-self.parameters['C'], self.parameters['C']] ])

class L(Circuit):
    """Inductor"""
    def __init__(self, plus, minus, L=0.0):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.branches.append(Branch(plus, minus))
        self.parameters['L']=L
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return array([[0.0 , 0.0, 1.0],
                       [0.0 , 0.0, -1.0],
                       [1.0 , -1.0, 0.0]])
    def C(self, x):
        return array([[0.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0],
                       [0.0, 0.0, self.parameters['L']]])

class VS(Circuit):
    """Independent voltage source

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c.solvedc()
    array([[ 1.5   ],
           [-0.0015]])
    
    """
    def __init__(self, plus=None, minus=None, v=0.0):
        Circuit.__init__(self, terminals=['plus', 'minus'])
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.branches.append(Branch(plus, minus))
        self.parameters['v']=v
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return array([[0.0 , 0.0, 1.0],
                       [0.0 , 0.0, -1.0],
                       [1.0 , -1.0, 0.0]])

    def U(self):
        return array([[0.0, 0.0, -self.parameters['v']]]).T

class IS(Circuit):
    """Independent current source

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c.solvedc()
    array([[ 1.]])
    
    """
    def __init__(self, plus, minus, i=0.0):
        Circuit.__init__(self, terminals=['plus', 'minus'])
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['i']=i
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def U(self):
        return array([[self.parameters['i'], -self.parameters['i']]]).T

class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2.0)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> c.solvedc()
    array([[ 1.5  ],
           [ 3.   ],
           [ 0.   ],
           [ 0.003]])
    """
    def __init__(self, inp, inn, outp, outn, g=1.0):
        Circuit.__init__(self, terminals=['inp','inn','outp','outn'])

        self.connectNode('inp', inp)
        self.connectNode('inn', inn)
        self.connectNode('outp', outp)
        self.connectNode('outn', outn)

        self.branches.append(Branch(outp, outn))

        self.parameters['g']=g
        self.terminals += (Terminal("inp"), Terminal("inn"), Terminal("outp"), Terminal("outn"))
        
    def G(self, x):
        return array([[0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 0.0, 0.0,-1.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0],
                      [self.parameters['g']  ,-self.parameters['g']  , -1.0,1.0, 0.0]])

class VCCS(Circuit):
    """Voltage controlled current source

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vccs'] = VCCS(n1, gnd, n2, gnd, gm=1e-3)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> c.solvedc()
    array([[ 1.5],
           [-1.5],
           [ 0. ]])

    """
    def __init__(self, inp, inn, outp, outn, gm=1e-3):
        Circuit.__init__(self, terminals=['inp','inn','outp','outn'])

        self.connectNode('inp', inp)
        self.connectNode('inn', inn)
        self.connectNode('outp', outp)
        self.connectNode('outn', outn)

        self.parameters['gm']=gm
        self.terminals += (Terminal("inp"), Terminal("inn"), Terminal("outp"), Terminal("outn"))
        
    def G(self, x):
        gm=self.parameters['gm']
        return array([[0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0],
                      [gm , -gm, 0.0, 0.0],
                      [-gm,  gm, 0.0, 0.0]])



gnd = Node("gnd")

def removeRowCol(A, n):
    for axis in range(len(A.shape)):
        A=delete(A, [n], axis=axis)
    return A

if __name__ == "__main__":
    import doctest
    doctest.testmod()
