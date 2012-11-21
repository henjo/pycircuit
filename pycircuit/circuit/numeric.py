# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric  operations that can be used as a toolkit for
Analysis objects

The module is based on `numpy <http://numpy.org>`_.

"""

from constants import *

import scipy.sparse as sp
import numpy as np
from numpy import cos, sin, tan, cosh, sinh, tanh, log, exp, pi, linalg,\
     inf, ceil, floor, linspace, eye, concatenate, sqrt, real, imag,\
     ones, complex, diff, delete, alltrue, maximum, size, conj

symbolic = False

ac_u_dtype = np.complex

def linearsolver(*args, **kvargs):
    return np.linalg.solve(*args, **kvargs)

def linearsolverError(*args, **kvargs):
    return np.linalg.LinAlgError

def toMatrix(array): 
    return array.astype('complex')

def det(x): 
    return np.linalg.det(x)

def dot(a,b):
    return a.dot(b)

def simplify(x): return x

def zeros(*args, **kvargs): 
    shape = args[0]
    
    ## Use numpy array vectors and sparse matrices
#    if type(shape) == int or len(shape) < 2:
    return np.zeros(*args, **kvargs)
#    else:
#        return sp.csr_matrix(shape, **kvargs)

def array(*args, **kvargs): 
    arr = np.array(*args, **kvargs)

    ## Use numpy array vectors and sparse matrices
    if len(arr.shape) < 2:
        return arr
    else:
        return sp.coo_matrix(arr, **kvargs)

def todense(a):
    return a.todense()

def array(*args, **kvargs):
    return np.array(*args, **kvargs)

def inv(*args, **kvargs): 
    return np.linalg.inv(*args, **kvargs)

def integer(x):
    return int(x)

numeric = True
    
