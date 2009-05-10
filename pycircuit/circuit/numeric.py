# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric  operations that can be used as a toolkit for
Analysis objects

The module is based on `numpy <http://numpy.org>`_.

"""

import numpy as np
from numpy import cos, sin, tan, cosh, sinh, tanh, log, exp, pi, linalg,\
     inf, ceil, floor, dot, linspace, eye, concatenate, sqrt, real, imag

from numpy.linalg import inv

ac_u_dtype = np.complex

def linearsolver(*args):
    args = [A.astype('complex') for A in args]

    return linalg.solve(*args)

def toMatrix(array): 
    return array.astype('complex')

def det(x): 
    return linalg.det(x)

def simplify(x): return x

def zeros(shape, dtype=None): return np.zeros(shape, dtype=dtype)
def array(x, dtype=None): return np.array(x, dtype=dtype)

numeric = True
    
