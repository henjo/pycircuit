# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric  operations that can be used as a toolkit for
Analysis objects

The module is based on `numpy <http://numpy.org>`_.

"""

from .constants import *

import numpy as np
from numpy import cos, sin, tan, cosh, sinh, tanh, log, exp, pi, linalg,\
     inf, ceil, floor, dot, linspace, eye, concatenate, sqrt, real, imag,\
     ones, diff, delete, alltrue, maximum, size, conj, cdouble

symbolic = False

ac_u_dtype = np.cdouble

def linearsolver(*args, **kvargs):
    return np.linalg.solve(*args, **kvargs)

def linearsolverError(*args, **kvargs):
    return np.linalg.LinAlgError

def toMatrix(array): 
    return array.astype('cdouble')

def det(x): 
    return np.linalg.det(x)

def simplify(x): return x

def zeros(*args, **kvargs): 
    return np.zeros(*args, **kvargs)

def array(*args, **kvargs): 
    return np.array(*args, **kvargs)

def inv(*args, **kvargs): 
    return np.linalg.inv(*args, **kvargs)

def integer(x):
    return int(x)

def complex(x):
    return complex(x)

numeric = True
    
