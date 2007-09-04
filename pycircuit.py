import unittest
from numpy import *
import scipy

class Branch:
    def __init__(self, plus, minus):
        self.plus = pl
        self.minus = node2

class Node:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name

class Terminal:
    def __init__(self, name):
        self.name = name

class Circuit:
    """The circuit class describes a full circuit, subcircuit or a single component. It contains
       a list of nodes,terminals and parameters.
       The terminals joins nodes inside the circuit with the nodes outside the circuit. When the circuit is
       instanciated, the outside nodes are passed to the object via the terminals.
       
       Attributes:
         nodes      A list that contains Node objects. If the circuit does not contain any internal nodes
                    the length is the same as the number of terminals.
         branches   A list of Branch objects. The solver will solve for the currents through the branches.
         terminals  A list that contains Terminal objects
         parameters A dictionary of parameter values keyed by the parameter name

    """
    def __init__(self):
        self.nodes=[]
        self.terminals=[]
        self.branches=[]
        self.parameters={}
        
    def G(self, x):
        """Calculate the G ((trans)conductance) matrix of the circuit given the x-vector"""
        pass

    def C(self, x):
        """Calculate the C ((trans)capacitance) matrix of the circuit given the x-vector"""
        pass

class SubCircuit(Circuit):
    """
    SubCircuit is aontainer for circuits.
    Attributes:
      elements        list of Circuit objects
      elementnodemap  list of translations between node indices of the elements to the
                      node index in the SubCircuit object.
    """
    def __init__(self):
        Circuit.__init__(self)
        self.elements = []
        self.elementnodemap = []
        
    def append(self, element):
        self.elements.append(element)

        # Add nodes from new element
        for node in element.nodes:
            if not node in self.nodes:
                self.nodes.append(node)
        print self.nodes
        
    def addNode(self, name=None):
        """Create an internal node in the circuit and return
           the new node."""
        newnode = Node(name)
        self.nodes.append(newnode)
        return newnode

    def updateNodeMap(self):
        """Update the elementnodemap attribute"""

        self.elementnodemap = []
        for element in self.elements:
            self.elementnodemap.append([self.nodes.index(node) for node in element.nodes])

    def G(self, x):
        N=len(self.nodes)+len(self.branches)
        G=zeros((N,N))

        self.updateNodeMap()
        
        for element, nodemap in zip(self.elements, self.elementnodemap):
            G[[[i] for i in nodemap], nodemap] += element.G(x)

        return G

class R(Circuit):
    """Resistor"""
    def __init__(self, plus, minus, R=1e3):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.parameters['R']=R
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return matrix([[1/self.parameters['R'] -1/self.parameters['R']],[-1/self.parameters['R'], 1/self.parameters['R']]])
        

class SimpleTests(unittest.TestCase):
    def testGMatrix(self):
        cir=SubCircuit()
        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")

        res = 1e3
        cir.append(R(net1, net2, res))

        G = cir.G(array([0,0]))
        self.assertEqual(G[0,0], 1/res)
        self.assertEqual(G[0,1], -1/res)
        print G
       
if __name__ == "__main__":
   unittest.main()

