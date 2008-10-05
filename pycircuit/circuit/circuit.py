# -*- coding: latin-1 -*-

from numpy import array, delete, linalg, size, zeros, concatenate, pi, dot, exp
from pycircuit.param import Parameter, ParameterDict
from constants import *
from copy import copy

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

    def __repr__(self):
        return 'Branch('+str(self.plus)+','+str(self.minus)+')'

### Default reference node
gnd = Node("gnd")

defaultepar = ParameterDict(
    Parameter("T", "Temperature", unit="K", default = 300.0))

class Circuit(object):
    """The circuit class describes a full circuit, subcircuit or a single component. 

    It contains a list of nodes,terminals, model and instance parameters.
    The terminals connect nodes inside the circuit to the nodes outside the circuit. When the circuit is
    instanciated, the outside nodes are passed to the object via the terminals.
       
    Attributes:
    nodes      -- A list that contains Node objects. If the circuit does not contain any internal nodes
                  the length is the same as the number of terminals.
    branches   -- A list of Branch objects. The solver will solve for the currents through the branches.
    terminals  -- A list that contains terminal names
    instparams -- A list of valid instance parameters (Parameter objects)
    mpar       -- A class variable with a ParameterDict containing model specific parameters
    ipar       -- A ParameterDict containing instance specific parameters
    nodenames  -- A dictionary that maps a local node name to the node object in nodes. If the node is
                  connnected to superior hierarchy levels through a terminal the terminal name must
                  be the same as the local node name

    """
    terminals = []
    mpar = ParameterDict()
    instparams = []
    def __init__(self, *args, **kvargs):
        self.nodes = []
        self.nodenames = {}
        self.branches = []
        self.ipar = ParameterDict(*self.instparams, **kvargs)
        self.x = {}

        for terminal in self.terminals:
            self.addNode(terminal)

        self.connectTerminals(**dict(zip(self.terminals, args)))
        
    def __copy__(self):
        newc = self.__class__()
        newc.nodes = copy(self.nodes)    
        newc.nodenames = copy(self.nodenames)    
        newc.branches = copy(self.branches)    
        newc.ipar = copy(self.ipar)
        newc.x = copy(self.x)
        return newc

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

    def getTerminalBranch(self, terminalname):
        """Find the branch that is connected to the given terminal of the given element.

        If no branch is found if there are more branches than one, None is returned

        Returns
        -------
        Tuple of Branch object and an integer indicating if terminal is connected to the positive or negative side
        of the branch. 1 == positive side, -1 == negative side
        
        >>> net1 = Node('net1')
        >>> net2 = Node('net2')
        >>> VS = VS(net1, net2)
        >>> VS.getTerminalBranch("minus")
        (Branch(Node('net1'),Node('net2')), -1)
        
        """
        plusbranches = [] ## Branches with its positive side connected to the terminal
        minusbranches = [] ## Branches with its negative side connected to the terminal
        for branch in self.branches:
            if branch.plus == self.nodenames[terminalname]:
                plusbranches.append(branch)
            elif branch.minus == self.nodenames[terminalname]:
                minusbranches.append(branch)

        if len(plusbranches + minusbranches) != 1:
            return None
        elif len(plusbranches) == 1:
            return plusbranches[0], 1
        elif len(minusbranches) == 1:
            return minusbranches[0], -1            

    def getNodeIndex(self, node):
        """Get row in the x vector of a node voltage"""
        return self.nodes.index(node)

    def getBranchIndex(self, branch):
        """Get row in the x vector of a branch current"""        
        
        return len(self.nodes) + self.branches.index(branch)

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
        
    def connectTerminals(self, **kvargs):
        """Connect nodes to terminals by using keyword arguments

        """
        for terminal, node in kvargs.items():
            if not terminal in self.terminals:
                raise ValueError('terminal '+str(terminal)+' is not defined')
            if node != None:
                self.nodes.remove(self.nodenames[terminal])
                if not node in self.nodes:
                    self.nodes.append(node)
                self.nodenames[terminal] = node

    @property
    def n(self):
        """Return size of x vector"""
        return len(self.nodes) + len(self.branches)

    def G(self, x, epar=defaultepar):
        """Calculate the G ((trans)conductance) matrix of the circuit given the x-vector"""
        return zeros((self.n, self.n), dtype=object)

    def C(self, x, epar=defaultepar):
        """Calculate the C (transcapacitance) matrix of the circuit given the x-vector"""
        return zeros((self.n, self.n), dtype=object)

    def U(self, t=0.0, epar=defaultepar):
        """Calculate the U column-vector the circuit at time t for a given x-vector.

        """
        return zeros((self.n, 1), dtype=object)

    def i(self, x, epar=defaultepar):
        """Calculate the i vector as a function of the x-vector

        For linear circuits i(x(t)) = G*x
        """
        return dot(self.G(x), x)

    def q(self, x, epar=defaultepar):
        """Calculate the q vector as a function of the x-vector

        For linear circuits q(x(t)) = C*x
        """
        return dot(self.G(x), x)

    def CY(self, x, epar=defaultepar):
        """Calculate the noise sources correlation matrix

        @param x:  the state vector
        @type  x:  numpy array
        @param epar: Environment parameters
        @type  epar: ParameterDict
        """
        return zeros((self.n, 1), dtype=object)

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

    def __repr__(self):
        return self.__class__.__name__ + \
               '(' + \
               ','.join([str(self.nodenames[term]) for term in self.terminals] + \
                        ['%s=%s'%(par.name, self.ipar.get(par)) for par in self.ipar.parameters]) + \
                        ')'
        
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
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.elements = {}
        self.elementnodemap = []
        self.elementbranchmap = []
        
    def __copy__(self):
        newc = super(SubCircuit, self).__copy__()
        newc.elements = copy(self.elements)
        newc.elementnodemap = copy(self.elementnodemap)
        newc.elementbranchmap = copy(self.elementbranchmap)
        return newc
    
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

    def getTerminalBranch(self, terminalname):
        """Find the branch that is connected to the given terminal of the given element.

        If no branch is found if there are more branches than one, None is returned

        The name can be a hierachical name and the notation is 'I1.I2.plus' for
        the terminal 'plus' of instance I2 of instance I1

        Returns
        -------
        Tuple of Branch object and an integer indicating if terminal is connected to the positive or negative side
        of the branch. 1 == positive side, -1 == negative side

        >>> c = SubCircuit()
        >>> net1 = c.addNode('net1')
        >>> net2 = c.addNode('net2')
        >>> c['vs'] = VS(net1, net2)
        >>> c.getTerminalBranch("vs.minus")
        (Branch(Node('net1'),Node('net2')), -1)
        
        """
        hierlevels = [part for part in terminalname.split('.')]

        if len(hierlevels)==1:
            Circuit.getTerminalBranch(self, hierlevels[0])
        else:
            return self.elements[hierlevels[0]].getTerminalBranch('.'.join(hierlevels[1:]))

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

    def G(self, x, epar=defaultepar):
        return self._add_element_submatrices('G', x, (epar,))

    def C(self, x, epar=defaultepar):
        return self._add_element_submatrices('C', x, (epar,))

    def U(self, t=0.0, epar=defaultepar):
        return self._add_element_subvectors('U', None, (t,epar))

    def i(self, x, epar=defaultepar):
        return self._add_element_subvectors('i', x, (epar,))

    def CY(self, x, epar=defaultepar):
        """Calculate composite noise source correlation matrix

        The noise sources in one element are assumed to be uncorrelated with the noise sources in the other elements.

        """
        return self._add_element_submatrices('CY', x, (epar,))
        
    def _add_element_submatrices(self, methodname, x, args):
        n=self.n
        A=zeros((n,n), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            if x != None:
                subx = x[nodemap,:]
                try:
                    A[[[i] for i in nodemap], nodemap] += getattr(element, methodname)(subx, *args)
                except Exception, e:
                    raise e.__class__(str(e) + ' at element '+str(element)+', args='+str(args))
            else:
                A[[[i] for i in nodemap], nodemap] += getattr(element, methodname)(*args)
        return A

    def _add_element_subvectors(self, methodname, x, args):
        n=self.n
        A=zeros((n,1), dtype=object)

        for instance,element in self.elements.items():
            nodemap = self.elementnodemap[instance]
            if x != None:
                subx = x[nodemap,:]
                A[nodemap,:] += getattr(element, methodname)(subx, *args)
            else:
                A[nodemap,:] += getattr(element, methodname)(*args)
        return A
        
        
class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['R']
    R(Node('1'),Node('gnd'),r=1000.0)
    >>> c.G(zeros((2,1)))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)
    >>> c = SubCircuit()
    >>> n2=c.addNode('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(zeros((2,1)))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm', default=1e3)]

    def G(self, x, epar=defaultepar):
        g = 1/self.ipar.r
        return  array([[g, -g],
                        [-g, g]], dtype=object)

    def CY(self, x, epar=defaultepar):
        if 'kT' in epar:
            iPSD = 4*epar.kT/self.ipar.r
        else:
            iPSD = 4*kboltzmann*epar.T/self.ipar.r
        return  array([[iPSD, -iPSD],
                       [-iPSD, iPSD]], dtype=object)
        

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(zeros((2,1)))
    array([[0, 0],
           [0, 0]], dtype=object)
    >>> c.C(zeros((2,1)))
    array([[1e-12, -1e-12],
           [-1e-12, 1e-12]], dtype=object)

    """

    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='c', desc='Capacitance', unit='F', default=1e-12)]

    def C(self, x, epar=defaultepar):
        return array([[self.ipar.c, -self.ipar.c],
                      [-self.ipar.c, self.ipar.c]])

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
    >>> c.G(zeros((3,1)))
    array([[0.0, 0.0, 1.0],
           [0.0, 0.0, -1.0],
           [1.0, -1.0, 0.0]], dtype=object)
    >>> c.C(zeros((3,1)))
    array([[0, 0, 0],
           [0, 0, 0],
           [0, 0, 1e-09]], dtype=object)
    """
    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='L', desc='Inductance', unit='H', default=1e-9)]

    _G = array([[0.0 , 0.0, 1.0],
                [0.0 , 0.0, -1.0],
                [1.0 , -1.0, 0.0]])
    def __init__(self, plus, minus, L=0.0):
        Circuit.__init__(self, plus, minus, L=L)
        self.branches.append(Branch(plus, minus))
    def G(self, x, epar):
        return self._G
    def C(self, x, epar):
        n = self.n
        C = zeros((n,n), dtype=object)
        C[-1,-1] = self.ipar.L
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
    instparams = [Parameter(name='v', desc='Source voltage', unit='V', default=1.0),
                  Parameter(name='noisePSD', desc='Voltage noise power spectral density', unit='V^2/Hz', default=0.0)]

    def __init__(self, plus, minus, **kvargs):
        Circuit.__init__(self, plus, minus, **kvargs)
        self.branches.append(Branch(plus, minus))

    def G(self, x, epar=defaultepar):
        return array([[0.0 , 0.0, 1.0],
                       [0.0 , 0.0, -1.0],
                       [1.0 , -1.0, 0.0]], dtype=object)

    def U(self, t=0.0, epar=defaultepar):
        return array([[0.0, 0.0, -self.ipar.v]], dtype=object).T

    def CY(self, x, epar=defaultepar):
        CY = super(VS, self).CY(x)
        CY[2, 2] = self.ipar.noisePSD

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

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
    instparams = [Parameter(name='i', desc='DC Current', unit='A', default=1e-3),
                  Parameter(name='noisePSD', desc='Current noise power spectral density', unit='A^2/Hz', default=0.0)]
    terminals = ['plus', 'minus']

    def U(self, t=0.0, epar=defaultepar):
        return array([[self.ipar.i, -self.ipar.i]]).T

    def CY(self, x, epar=defaultepar):
        return  array([[self.ipar.noisePSD, -self.ipar.noisePSD],
                       [-self.ipar.noisePSD, self.ipar.noisePSD]], dtype=object)

class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2.0)
    >>> c['vcvs'].nodes
    [Node('2'), Node('gnd'), Node('1')]
    >>> c['vcvs'].branches
    [Branch(Node('2'),Node('gnd'))]
    >>> c['vcvs'].G(zeros((4,1)))
    array([[0, 0, 0, 1.0],
           [0, 0, 0, -1.0],
           [0, 0, 0, 0],
           [-1.0, -1.0, 2.0, 0]], dtype=object)
    """
    instparams = [Parameter(name='g', desc='Voltage gain', unit='V/V', default=1.0)]
    terminals = ('inp', 'inn', 'outp', 'outn')
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], self.nodenames['outn']))
        
    def G(self, x, epar=defaultepar):
        G = super(VCVS, self).G(x)
        branchindex = -1
        inpindex,innindex,outpindex,outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, branchindex] += 1.0
        G[outnindex, branchindex] += -1.0
        G[branchindex, outpindex] += -1.0
        G[branchindex, outnindex] += 1.0
        G[branchindex, inpindex] += self.ipar.g
        G[branchindex, innindex] += -self.ipar.g
        return G

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
    instparams = [Parameter(name='gm', desc='Transconductance', unit='A/V', default=1e-3)]
    
    def G(self, x, epar=defaultepar):
        G = super(VCCS, self).G(x)
        gm=self.ipar.gm
        inpindex,innindex,outpindex,outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, inpindex] += gm
        G[outpindex, innindex] += -gm
        G[outnindex, inpindex] += -gm
        G[outnindex, innindex] += gm
        return G

class Nullor(Circuit):
    """Nullor

    From Wikipedia, the free encyclopedia

     A nullor is a theoretical two-port network comprised of a nullator at its input and a norator at its output.[1]
     Nullors represent an ideal amplifier, having infinite current, voltage, transconductance and transimpedance gain.[2] Its
     transmission parameters are all zero.

      
     1. The name "nullor" was introduced by H.J. Carlin, Singular network elements, IEEE Trans. Circuit Theory, March 1965, vol. CT-11, pp. 67-72.
 
     2. Verhoeven C J M van Staveren A Monna G L E Kouwenhoven M H L & Yildiz E (2003). Structured electronic design: negative feedback amplifiers.
          Boston/Dordrecht/London: Kluwer Academic, §2.2.2 pp. 32-34. ISBN 1402075901.

    """
    terminals = ('inp', 'inn', 'outp', 'outn')

    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Nullor, self).G(x)
        branchindex = -1
        inpindex,innindex,outpindex,outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in ('inp', 'inn', 'outp', 'outn'))

        G[outpindex, branchindex] += 1.0   
        G[outnindex, branchindex] += -1.0
        G[branchindex, inpindex] += 1.0
        G[branchindex, innindex] += -1.0
        return G


class Transformer(Circuit):
    """Ideal transformer

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.addNode('1')
    >>> n2=c.addNode('2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = Transformer(n1, gnd, n2, gnd, n=2.0)
    >>> c['vcvs'].nodes
    [Node('2'), Node('gnd'), Node('1')]
    >>> c['vcvs'].branches
    [Branch(Node('2'),Node('gnd'))]
    >>> c['vcvs'].G(zeros((4,1)))
    array([[0, 0, 0, 1.0],
           [0, 0, 0, -3.0],
           [0, 0, 0, 2.0],
           [2.0, -1.0, -1.0, 0]], dtype=object)
    """
    instparams = [Parameter(name='n', desc='Winding ratio', unit='', default=1.0)]
    terminals = ('inp', 'inn', 'outp', 'outn')
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Transformer, self).G(x)
        branchindex = -1
        inpindex,innindex,outpindex,outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in ('inp', 'inn', 'outp', 'outn'))
        G[inpindex, branchindex] += self.ipar.n
        G[innindex, branchindex] += -self.ipar.n
        G[outpindex, branchindex] += 1.0
        G[outnindex, branchindex] += -1.0
        G[branchindex, outpindex] += self.ipar.n
        G[branchindex, outnindex] += -self.ipar.n
        G[branchindex, inpindex] += -1.0
        G[branchindex, innindex] += 1.0
        return G

class Diode(Circuit):
    terminals = ['plus', 'minus']
    mpar = Circuit.mpar.copy( Parameter(name='IS', desc='Saturation current', unit='A', default=1e-13) )
        
    def G(self, x, epar=defaultepar):
        VD = x[0,0]-x[1,0]
        VT = kboltzmann*epar.T / qelectron
        g = self.mpar.IS*exp(VD/VT)/VT
        return array([[g, -g],
                      [-g, g]], dtype=object)

    def i(self, x, epar=defaultepar):
        """
        
        """
        VD = x[0,0]-x[1,0]
        VT = kboltzmann*epar.T / qelectron
        I = self.mpar.IS*(exp(VD/VT)-1.0)
        return array([[I, -I]], dtype=object).T

if __name__ == "__main__":
    import doctest
    doctest.testmod()
