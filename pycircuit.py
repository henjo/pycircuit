import unittest

class Node:
    def __init__(self, name):
        self.name = name

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
         terminals  A list that contains Terminal objects
         parameters A dictionary of parameter values keyed by the parameter name

    """
    def __init__(self):
        self.nodes=[]
        self.terminals=[]
        self.parameters={}
        
    def G(self, x):
        """Calculate the G ((trans)conductance) matrix of the circuit given the x-vector"""
        pass

    def C(self, x):
        """Calculate the C ((trans)capacitance) matrix of the circuit given the x-vector"""
        pass

class SubCircuit(Circuit):
    def __init__(self):
        Circuit.__init__(self)
        self.elements = []
        
    def append(self, circuit):
        self.elements.append(circuit)

    def addNode(self, name=None):
        """Create an internal node in the circuit and return
           the new node."""
        return self.nodes.append(Node(name))

class Resistor(Circuit):
    def __init__(self, plusnet, minusnet, R=1e3):
        Circuit.__init__(self)
        self.parameters['R']=R
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
class SimpleTests(unittest.TestCase):
    def testGMatrix(self):
        cir=SubCircuit()
        net1 = cir.addNode()
        net2 = cir.addNode()

        cir.append(Resistor(net1, net2, 1e3))
       
if __name__ == "__main__":
   unittest.main()

