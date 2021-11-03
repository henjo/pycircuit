# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

#import yapgvb
import os
from copy import copy
import sympy

class DDD(object):
    pass

class Node(object):
    def __init__(self, index, D1=None, D0=None, sign=1):
        """Construct a DDD (Determinant-Decision-Diagram)
        from top vertex values and 1- and 0-edges

        @param index: vertex index
        @param D1:  1-edge 
        @type  D1:  Node instance
        @param D0:  0-edge 
        @type  D1:  Node instance
        @param sign: sign of vertex (+1, -1)
        @type  sign: int
        
        """
        self.index = index
        self.sign = sign

        # Check if terminal nodes
        if index in (0,1):
            self.D1 = None
            self.D0 = None
        else:
            # Default arguments
            if D1 is None:
                self.D1 = VertexOne
            else:
                self.D1 = D1

            if D0 is None:
                self.D0 = VertexZero
            else:
                self.D0 = D0
        
    def cofactor(self, s):
        """Return cofactor of DDD with regards to s"""
        if self.index < s:
            return Node(0)
        elif self.index == s:
            return self.D1
        elif self.index > s:
            return Node(self.index, D1=self.D1.cofactor(s), D0=self.D0.cofactor(s))

    def remainder(self, s):
        if self.index < s:
            return self
        elif self.index == s:
            return self.D0
        elif self.index > s:
            return Node(self.index, D1=self.D1.remainder(s), D0=self.D0.remainder(s))

    def isleaf(self):
        return self.D0 is None and self.D1 is None

    def __eq__(self, P):
        if not isinstance(P, Node):
            return False
        else:
            return self.index == P.index and self.D0 == P.D0 and self.D1 == P.D1 \
                   and self.sign == P.sign
        
    def union(self, P):
        if self == VertexZero:
            return P
        elif P == VertexZero:
            return self
        elif self == P:
            return self
        elif self.index > P.index:
            return Node(self.index, D1=self.D1, D0=self.D0.union(P))
        elif self.index < P.index:
            return Node(P.index, D1=P.D1, D0=self.union(P.D0))
        elif self.index == P.index:
            return Node(self.index, D1=self.D1.union(P.D1), D0=self.D0.union(P.D0))

    def intersec(self, P):
        if self == VertexZero or P == VertexZero:
            return VertexZero
        elif self == P:
            return self
        elif self.index > P.index:
            return self.D0.intersec(P)
        elif self.index < P.index:
            return self.intersec(P.D0)
        elif self.index == P.index:
            return Node(self.index, D1=self.D1.intersec(P.D1), D0=self.D0.intersec(P.D0))

    def __mul__(self, s):
        if self.index < s:
            return Node(s, D1=self, D0=VertexZero)
        elif self.index == s:
            return copy(self)
        elif self.index > s:
            return Node(self.index, D1=self.D1*s, D0=self.D0*s)

    def __sub__(self, P):
        if self.index == 0 and self.isleaf():
            return VertexZero
        elif P.index == 0 and P.isleaf():
            return self
        elif self == P:
            return VertexZero
        elif self.index > P.index:
            return Node(self.index, D1=self.D1, D0=self.D0-P)
        elif self.index < P.index:
            return self - P.D0
        elif self.index == P.index:
            return Node(self.index, D1=self.D1 - P.D1, D0=self.D0-P.D0)
        
    def __or__(self, s):
        return self.union(s)

    def __and__(self, s):
        return self.intersec(s)

    def eval(self):
        if self.D1 is None and self.D0 is None:
            return self.index
        else:
            return self.D0.eval() + self.sign * self.index * self.D1.eval()

    def __repr__(self):
        return str(self.eval())

    def _asdotAddSelf(self, g):
        if self.isleaf():
            top = g.add_node(str(self), label=str(self.index), shape='box')
        else:
            top = g.add_node(str(self), label=str(self.index))

            D0top = self.D0._asdotAddSelf(g)
            D1top = self.D1._asdotAddSelf(g)

            top >> D1top
            D1edge = top >> D0top
            D1edge.style='dashed'
            
        return top
        
    def asdot(self, name=None):
        """Return GraphViz representation of Node"""
        g = yapgvb.Digraph()
        topnode = self._asdotAddSelf(g)

        if name:
            namenode = g.add_node(label=name, shape='none')
            namenode >> topnode
            
        return g

VertexOne = Node(1, None, None)
VertexZero = Node(0, None, None)

def DDD_of_matrix(A):
    """Create DDD of a Sympy matrix
    @param A:  Input matrix
    @type A:   sympy.Matrix

    >>> a,b,c,d,e,f,g,h,i,j=sympy.symbols('abcdefghij')
    >>> A = sympy.Matrix((a,b,0,0),(c,d,e,0),(0,f,g,h),(0,0,i,j))
    >>> DDD_matrix(A)
    
    """
    print A
    if A == sympy.zero(1):
        return VertexZero
    if A == sympy.one(1):
        return VertexOne

    # Select row, column
    i = 0; j = 0
    D0 = copy(A)
    D0[i,j] = 0
    return DDD(A[i,j],
               D1 = DDD_of_matrix(A.minorMatrix(i,j)),
               D0 = D0)    

if __name__ == '__main__':
    from sympy import Symbol
    
    A = VertexZero
    B = VertexOne
    C = B * Symbol('x1')
    D = B * Symbol('x2')
    E = C * Symbol('x2')
    F = C | E
    G = E.cofactor(Symbol('x2'))
    H = F.remainder(Symbol('x2'))
    I = F - C
    J = F & C
    
    print ','.join(map(str, [A,B,C,D,E,F,G,H,I,J]))

#    def displaygraph(D):
#        D.asdot().write('/tmp/graph.dot')
#        os.system('dot -Txlib /tmp/graph.dot')

#    displaygraph(Eg)
    
 
