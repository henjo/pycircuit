# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric  operations that can be used as a toolkit for
Analysis objects

The module is based on `numpy <http://numpy.org>`_.

"""

from numpy import *

def linearsolver(*args):
    args = [A.astype('complex') for A in args]

    return linalg.solve(*args)

def toMatrix(array): 
    return array.astype('complex')

def det(x): 
    return linalg.det(x)

numeric = True
    
