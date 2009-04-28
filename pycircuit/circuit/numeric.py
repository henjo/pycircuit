# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric  operations that can be used as a toolkit for
Analysis objects

The module is based on `numpy <http://numpy.org>`_.

"""

import numpy
from numpy import cos, sin, tan, exp, pi, linalg, inf, ceil, floor, dot

def linearsolver(*args):
    args = [A.astype('complex') for A in args]

    return linalg.solve(*args)

def toMatrix(array): 
    return array.astype('complex')

def det(x): 
    return linalg.det(x)

def simplify(x): return x

def zeros(shape, dtype=float): return numpy.zeros(shape, dtype=dtype)
def array(x, dtype=float): return numpy.array(x, dtype=dtype)

numeric = True
    
