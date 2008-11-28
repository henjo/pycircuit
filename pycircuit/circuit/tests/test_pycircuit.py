# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from numpy import *
from scipy import *
from pycircuit import *
import unittest

class SimpleTests(unittest.TestCase):
    def testGMatrix(self):
        cir=SubCircuit()
        net1 = cir.add_node("net1")
        net2 = cir.add_node("net2")

        res = 1e3
        cir['R1'] = R(net1, net2, res)
        cir['R2'] = R(net1, net2, res)

        cir['VS'] = VS(net1, gnd, 1.0)

        G = cir.G(array([0,0]))
        U = cir.U()
        self.assertEqual(G[0,0], 2/res)
        self.assertEqual(G[0,1], -2/res)

    def testrc(self):
        cir=SubCircuit()
        net1 = cir.add_node("net1")
        net2 = cir.add_node("net2")

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
        net1 = cir.add_node("net1")
        net2 = cir.add_node("net2")

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

        net1 = cir.add_node("net1")
        net2 = cir.add_node("net2")
        
        v0,R1,R2=map(Symbol, ('v0','R1','R2'))

        cir['R1']=R(net1, net2, R1)
        cir['R2']=R(net2, gnd, R2)
        cir['VS']=VS(net1, gnd, v0)

        res = cir.solvesymbolic()
        self.assertEqual(simplify(res[1,0]-v0*R2/(R1+R2)), 0.0)

    def testRCfilter(self):
        
        cir=SubCircuit()

        net1 = cir.add_node("net1")
        net2 = cir.add_node("net2")
        
        v0,R1,C1=map(Symbol, ('v0','R1','C1'))

        cir['R1']=R(net1, net2, R1)
        cir['R2']=C(net2, gnd, C1)
        cir['VS']=VS(net1, gnd, v0)

        res = cir.solvesymbolic()
        self.assertEqual(res[1,0]-v0/(1+Symbol('s')*R1*C1), 0)

if __name__ == "__main__":
    unittest.main()
