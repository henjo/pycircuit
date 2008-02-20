import unittest
from numpy import *
import pylab
import sympy
from sympy import Symbol, Matrix, symbols, simplify

class Branch:
    def __init__(self, plus, minus, name=None):
        self.plus = plus
        self.minus = minus
        self.name = name

class Node:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name

class Terminal:
    def __init__(self, name):
        self.name = name

gnd = Node("gnd")

def removeRowCol(A, n):
    for axis in range(len(A.shape)):
        A=delete(A, [n], axis=axis)
    return A

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

    def addNode(self, name=None):
        """Create an internal node in the circuit and return
           the new node."""
        newnode = Node(name)
        self.nodes.append(newnode)
        return newnode

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
    """Resistor"""
    def __init__(self, plus, minus, R=1e3):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.parameters['R']=R
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return array([[1/self.parameters['R'] , -1/self.parameters['R']],
                       [-1/self.parameters['R'], 1/self.parameters['R']] ])

class C(Circuit):
    """Capacitor"""
    def __init__(self, plus, minus, C=0.0):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.parameters['C']=C
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def C(self, x):
        return array([[self.parameters['C'] , -self.parameters['C']],
                       [-self.parameters['C'], self.parameters['C']] ])

class VS(Circuit):
    """Independent voltage source"""
    def __init__(self, plus, minus, V=0.0):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.branches.append(Branch(plus, minus))
        self.parameters['V']=V
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return array([[0.0 , 0.0, 1.0],
                       [0.0 , 0.0, -1.0],
                       [1.0 , -1.0, 0.0]])

    def U(self):
        return array([[0.0, 0.0, -self.parameters['V']]]).T

class VCVS(Circuit):
    """Voltage controlled voltage source"""
    def __init__(self, inp, inn, outp, outn, g=1.0):
        Circuit.__init__(self)
        map(self.nodes.append, [inp,inn,outp,outn])
        self.branches.append(Branch(outp, outn))
        self.parameters['g']=g
        self.terminals += (Terminal("inp"), Terminal("inn"), Terminal("outp"), Terminal("outn"))
        
    def G(self, x):
        return array([[0.0, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 0.0, 0.0,-1.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0],
                      [g  ,-g  , 1.0,-1.0, 0.0]])

class VCCS(Circuit):
    """Voltage controlled current source"""
    def __init__(self, inp, inn, outp, outn, gm=1e-3):
        Circuit.__init__(self)
        map(self.nodes.append, [inp,inn,outp,outn])
        self.parameters['gm']=gm
        self.terminals += (Terminal("inp"), Terminal("inn"), Terminal("outp"), Terminal("outn"))
        
    def G(self, x):
        return array([[0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0],
                      [gm , -gm, 0.0, 0.0],
                      [-gm,  gm, 0.0, 0.0]])

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

class IS(Circuit):
    """Independent current source"""
    def __init__(self, plus, minus, I=0.0):
        Circuit.__init__(self)
        self.nodes.append(plus)
        self.nodes.append(minus)
        self.parameters['I']=I
        self.terminals += (Terminal("plus"), Terminal("minus"))
        
    def G(self, x):
        return array([[0.0 , 0.0],
                       [0.0 , 0.0]])

    def U(self):
        return array([self.parameters['I'], -self.parameters['I']]).T

class SimpleTests(unittest.TestCase):
    def testGMatrix(self):
        cir=SubCircuit()
        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")

        res = 1e3
        cir['R1'] = R(net1, net2, res)
        cir['R2'] = R(net1, net2, res)

        cir['VS'] = VS(net1, gnd, 1.0)

        G = cir.G(array([0,0]))
        U = cir.U()
        self.assertEqual(G[0,0], 2/res)
        self.assertEqual(G[0,1], -2/res)

        print cir.solvedc()

    def testrc(self):
        cir=SubCircuit()
        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")

        cir['R1'] = R(net1, net2, 1e3)
        cir['C1'] = C(net2, gnd, 1e-12)
        cir['VS'] = VS(net1, gnd, 1.0)

        f = logspace(6,9)
        result = array(cir.solveac(f))
#        pylab.semilogx(f, 20*log10(abs(result[:,1])))
#        pylab.grid()
#        pylab.show()

    def testlc(self):
        cir=SubCircuit()
        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")

        cir['L1'] = L(net1, net2, 1e-3)
        cir['C1'] = C(net2, gnd, 1e-12)
        cir['VS'] = VS(net1, gnd, 1.0)

        f = logspace(6,9)
        result = array(cir.solveac(f))
#        pylab.semilogx(f, 20*log10(abs(result[:,1])))
#        pylab.grid()
#        pylab.show()

class SymbolicTests(unittest.TestCase):
    def testvoltagedivider(self):
        
        cir=SubCircuit()

        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")
        
        v0,R1,R2=map(Symbol, ('v0','R1','R2'))

        cir['R1']=R(net1, net2, R1)
        cir['R2']=R(net2, gnd, R2)
        cir['VS']=VS(net1, gnd, v0)

        res = cir.solvesymbolic()
    
        self.assertEqual(simplify(res[Symbol('vnet2')]-v0*R2/(R1+R2)), 0.0)

    def testRCfilter(self):
        
        cir=SubCircuit()

        net1 = cir.addNode("net1")
        net2 = cir.addNode("net2")
        
        v0,R1,C1=map(Symbol, ('v0','R1','C1'))

        cir['R1']=R(net1, net2, R1)
        cir['R2']=C(net2, gnd, C1)
        cir['VS']=VS(net1, gnd, v0)

        res = cir.solvesymbolic()
        self.assertEqual(res[Symbol('vnet2')]-v0/(1+Symbol('s')*R1*C1), 0)
        print res

if __name__ == "__main__":
    unittest.main()

