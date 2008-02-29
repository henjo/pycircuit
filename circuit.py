import unittest
import numpy
from numpy import array, delete, linalg, size, zeros, concatenate, pi
import pylab
import sympy

class Node(object):
    """A node is a region in an electric circuit where the voltage is the same.
    
    A node can have a name
    """
    def __init__(self, name=None):
        self.name = name

    def __repr__(self):
        return self.__class__.__name__ + '(\'' + self.name + '\')'

class Branch(object):
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

### Default reference node
gnd = Node("gnd")

class Circuit(object):
    """The circuit class describes a full circuit, subcircuit or a single component. 

    It contains a list of nodes,terminals and parameters.
    The terminals connect nodes inside the circuit to the nodes outside the circuit. When the circuit is
    instanciated, the outside nodes are passed to the object via the terminals.
       
    Attributes:
    nodes      -- A list that contains Node objects. If the circuit does not contain any internal nodes
                  the length is the same as the number of terminals.
    branches   -- A list of Branch objects. The solver will solve for the currents through the branches.
    terminals  -- A list that contains terminal names
    parameters -- A dictionary of parameter values keyed by the parameter name
    nodenames  -- A dictionary that maps a local node name to the node object in nodes. If the node is
                  connnected to superior hierarchy levels through a terminal the terminal name must
                  be the same as the local node name
    


    """
    terminals = []
    def __init__(self):
        self.nodes = []
        self.nodenames = {}
        self.branches = []
        self.parameters = {}
        self.x = {}

        for terminal in self.terminals:
            self.addNode(terminal)
        
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

    def getNodeIndex(self, node):
        """Get row in the x vector of a node"""
        return self.nodes.index(node)

    def getNode(self, name):
        """Find a node by name.
        
        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.getNode('n1')
        Node('n1')
        
        """
        return self.nodenames[name]

    def getNodeName(self, node):
        """Find the name of a node
        
        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.getNodeName(n1)
        'n1'

        """
        
        for k, v in self.nodenames.items():
            if v == node:
                return k
        
    def connectNode(self, terminal, node):
        """Connect an external node to a terminal

        """
        if node != None:
            if not terminal in self.terminals:
                raise ValueError('terminal '+str(terminal)+' is not defined')
            self.nodes.remove(self.nodenames[terminal])
            self.nodes.append(node)
            self.nodenames[terminal] = node

    def n(self):
        """Return size of x vector"""
        return len(self.nodes) + len(self.branches)

    def G(self, x):
        """Calculate the G ((trans)conductance) matrix of the circuit given the x-vector"""
        return zeros((self.n(), self.n()), dtype=object)

    def C(self, x):
        """Calculate the C (transcapacitance) matrix of the circuit given the x-vector"""
        return zeros((self.n(), self.n()), dtype=object)

    def U(self, t=0.0):
        """Calculate the U (circuit input) column-vector the circuit at time t"""
        return zeros((self.n(), 1), dtype=object)

    def nameStateVector(self, x, analysis=''):
        """Map state variables with names and return the state variables in a dictionary keyed by the names

        >>> c = SubCircuit()
        >>> n1=c.addNode('net1')
        >>> c['is'] = IS(gnd, n1, i=1e-3)
        >>> c['R'] = R(n1, gnd, r=1e3)
        >>> c.nameStateVector(array([[1.0]]))
        {'net1': 1.0}

        >>> 

        """
        result = {}
        for xvalue, node in zip(x[:len(self.nodes)][0], self.nodes):
            result[self.getNodeName(node)] = xvalue

        for i, xvalue, branch in enumerate(zip(x[len(self.nodes):], self.branches)):
            result['i' + analysis + str(i) + ')'] = xvalue            

        return result

        
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

    def __getitem__(self, instancename):
        """Get instance"""
        return self.elements[instancename]

    def getNode(self, name):
        """Find a node by name.

        The name can be a hierachical name and the notation is 'I1.I2.net1' for
        a net net1 in instance I2 of instance I1
        
        >>> c = Circuit()
        >>> n1 = c.addNode("n1")
        >>> c.getNode('n1')
        Node('n1')

        >>> c1 = SubCircuit()
        >>> c2 = SubCircuit()
        >>> c1['I1'] = c2
        >>> n1 = c2.addNode("net1")
        >>> c1.getNode('I1.net1')
        Node('net1')
        
        """
        hierlevels = [part for part in name.split('.')]
            
        if len(hierlevels)==1:
            return self.nodenames[hierlevels[0]]
        else:
            return self.elements[hierlevels[0]].getNode('.'.join(hierlevels[1:]))

    def getNodeName(self, node):
        """Find the name of a node
        
        >>> c1 = SubCircuit()
        >>> c2 = SubCircuit()
        >>> c1['I1'] = c2
        >>> n1 = c2.addNode("net1")
        >>> c1.getNodeName(n1)
        'net1'

        """

        ## Use name of node object if present
        if node.name != None:
            return node.name
        
        ## First search among the local nodes
        name = Circuit.getNodeName(self, node)
        if name != None:
            return name
        
        ## Then search in the circuit elements
        for instname, element in self.elements.items():
            name =  element.getNodeName(node)
            if name != None:
                return instname + '.' + name
        
    def updateNodeMap(self):
        """Update the elementnodemap attribute"""

        self.elementnodemap = {}
        for instance, element in self.elements.items():
            self.elementnodemap[instance] = [self.nodes.index(node) for node in element.nodes] + \
                                            [self.branches.index(branch)+len(self.nodes) for branch in element.branches]

    def G(self, x):
        n=self.n()
        G=zeros((n,n), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            G[[[i] for i in nodemap], nodemap] += element.G(x)
            
        return G

    def C(self, x):
        n=self.n()
        C=zeros((n,n), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            C[[[i] for i in nodemap], nodemap] += element.C(x)
        return C

    def U(self, t=0.0):
        n=self.n()
        U=zeros((n,1), dtype=object)

        self.updateNodeMap()
        
        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            U[nodemap] += element.U(t)
        return U

class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c.G(0)
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)
    >>> c = SubCircuit()
    >>> n2=c.addNode('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(0)
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']

    def __init__(self, plus, minus, r=1e3):
        Circuit.__init__(self, )
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['r']=r
        
    def G(self, x):
        return array([[1/self.parameters['r'], -1/self.parameters['r']],
                      [-1/self.parameters['r'], 1/self.parameters['r']]], dtype=object)

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(0)
    array([[0, 0],
           [0, 0]], dtype=object)
    >>> c.C(0)
    array([[1e-12, -1e+12],
           [-1e+12, 1e+12]], dtype=object)

    """

    terminals = ['plus', 'minus']    
    def __init__(self, plus, minus, c=0.0):
        Circuit.__init__(self, )
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['c'] = c

    def C(self, x):
        return array([[self.parameters['c'], -1/self.parameters['c']],
                      [-1/self.parameters['c'], 1/self.parameters['c']]])

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['C'] = L(n1, gnd, L=1e-9)
    >>> c.G(0)
    array([[0.0, 0.0, 1.0],
           [0.0, 0.0, -1.0],
           [1.0, -1.0, 0.0]], dtype=object)
    >>> c.C(0)
    array([[0, 0, 0],
           [0, 0, 0],
           [0, 0, 1e-09]], dtype=object)
    """
    terminals = ['plus', 'minus']    

    _G = array([[0.0 , 0.0, 1.0],
                [0.0 , 0.0, -1.0],
                [1.0 , -1.0, 0.0]])
    def __init__(self, plus, minus, L=0.0):
        Circuit.__init__(self)
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.branches.append(Branch(plus, minus))
        self.parameters['L']=L
        
    def G(self, x):
        return self._G
    def C(self, x):
        n = self.n()
        C = zeros((n,n), dtype=object)
        C[-1,-1] = self.parameters['L']
        return C

class VS(Circuit):
    """Independent voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd)
    array([[ 1.5   ],
           [ 0.    ],
           [-0.0015]])
    
    """
    terminals = ['plus', 'minus']

    def __init__(self, plus=None, minus=None, v=0.0):
        Circuit.__init__(self)
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.branches.append(Branch(plus, minus))
        self.parameters['v']=v

    def G(self, x):
        return array([[0.0 , 0.0, 1.0],
                       [0.0 , 0.0, -1.0],
                       [1.0 , -1.0, 0.0]], dtype=object)

    def U(self, t=0.0):
        return array([[0.0, 0.0, -self.parameters['v']]], dtype=object).T

class IS(Circuit):
    """Independent current source

    >>> from analysis import DC, gnd as gnd2
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd)
    array([[ 1.],
           [ 0.]])
    
    """
    terminals = ['plus', 'minus']

    def __init__(self, plus, minus, i=0.0):
        Circuit.__init__(self)
        self.connectNode('plus', plus)
        self.connectNode('minus', minus)
        self.parameters['i']=i
        
    def U(self, t=0.0):
        return array([[self.parameters['i'], -self.parameters['i']]]).T

class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2.0)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd)
    array([[ 1.5  ],
           [ 3.   ],
           [ 0.   ],
           [ 0.   ],
           [ 0.003]])
    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    def __init__(self, inp, inn, outp, outn, g=1.0):
        Circuit.__init__(self)

        self.connectNode('inp', inp)
        self.connectNode('inn', inn)
        self.connectNode('outp', outp)
        self.connectNode('outn', outn)

        self.branches.append(Branch(outp, outn))

        self.parameters['g']=g
        
    def G(self, x):
        return array([[0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 0.0, 0.0,-1.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0],
                      [self.parameters['g']  ,-self.parameters['g']  , -1.0,1.0, 0.0]])

class VCCS(Circuit):
    """Voltage controlled current source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vccs'] = VCCS(n1, gnd, n2, gnd, gm=1e-3)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd)
    array([[ 1.5],
           [-1.5],
           [ 0. ],
           [ 0. ]])

    """
    terminals = ['inp', 'inn', 'outp', 'outn']
    
    def __init__(self, inp, inn, outp, outn, gm=1e-3):
        Circuit.__init__(self)

        self.connectNode('inp', inp)
        self.connectNode('inn', inn)
        self.connectNode('outp', outp)
        self.connectNode('outn', outn)

        self.parameters['gm']=gm
        
    def G(self, x):
        gm=self.parameters['gm']
        return array([[0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0],
                      [gm , -gm, 0.0, 0.0],
                      [-gm,  gm, 0.0, 0.0]])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
