# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.ddd import *
import unittest

class DDDBasicTests(unittest.TestCase):
    def testCreateInstance(self):
        """Test instantiation of DDD"""
        for D in (DDD('x1',D1=VertexOne, D0=VertexZero,sign=-1),
                  DDD('x1',D1=VertexOne, sign=-1),
                  DDD('x1',D0=VertexZero, sign=-1),
                  DDD('x1',sign=-1)
                  ):
            self.assertEqual(D.top, 'x1')
            self.assertEqual(D.D1.top, 1)
            self.assertEqual(D.D1.D1, None)
            self.assertEqual(D.D1.D0, None)
            self.assertEqual(D.D0.top, 0)
            self.assertEqual(D.D0.D1, None)
            self.assertEqual(D.D0.D0, None)
            self.assertEqual(D.sign, -1)
    def testEqualSimple(self):
        D = DDD(42,sign=-1)
        self.assertEqual(D, DDD(42, sign=-1))
        self.assertNotEqual(D, DDD(43,D1=VertexOne, D0=VertexZero,sign=-1))
        self.assertNotEqual(D, DDD(1,sign=-1))
        self.assertNotEqual(D, DDD(42,D1=VertexOne, D0=VertexZero,sign=1))
    def testEqualLessSimple(self):
        D1 = DDD(43)
        D = DDD(42,D1=D1, D0=VertexZero,sign=-1)
        self.assertEqual(D, DDD(42,D1=D1, D0=VertexZero,sign=-1))
        self.assertNotEqual(D, DDD(43,D1=VertexOne, D0=VertexZero,sign=-1))
        self.assertNotEqual(D, DDD(42,D1=VertexOne, D0=None,sign=-1))
        self.assertNotEqual(D, DDD(42,D1=None, D0=VertexZero,sign=-1))
        self.assertNotEqual(D, DDD(42,D1=VertexOne, D0=VertexZero,sign=1))
    def testEvaluate(self):
        D1 = DDD(42)
        D2 = DDD(42, sign=-1)
        self.assertEqual(D1.eval(), 42)
        self.assertEqual(D2.eval(), -42)

        D1 = DDD(43)
        D2 = DDD(42,D1=D1, D0=VertexZero,sign=-1)
        self.assertEqual(D2.eval(), -43*42)

class DDDTestFromBook(unittest.TestCase):
    """Test from "Symbolic Analysis and reduction of VLSI circuits" p.198
    """
    def test(self):
        A = VertexZero
        B = VertexOne
        C = B * 42
        D = B * 42
        E = C * 43
        F = C | E
        G = E.cofactor(43)
        H = F.remainder(43)
        I = F - C
        J = F & C

        self.assertEqual(C,G)
        self.assertEqual(C,H)
        self.assertEqual(C,J)
        self.assertEqual(E,I)

        
if __name__ == "__main__":
    unittest.main()
