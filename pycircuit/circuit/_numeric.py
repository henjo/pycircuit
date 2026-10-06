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
     ones, diff, delete, all, maximum, minimum, size, conj, cdouble, sum, max, where, abs, insert,\
     arctan2



def alltrue(a, *args, **kwargs):
    """numpy's `all` (imported above) -- by the array's own method where `a`
    is an exact ndarray and nothing else is asked: the same
    `logical_and.reduce(a, None, bool)` (`numpy.all` and `ndarray.all`
    reduce alike, result and type), without the module function's dispatch:
    4.5 k instructions against 11.6 k a call, 16.6 k calls a vdP PSS (the C
    lookup's and the Newton's convergence tests; speed round 12, stage 3)."""
    if type(a) is np.ndarray and not args and not kwargs:
        return a.all()
    return all(a, *args, **kwargs)

# natural logarithm; matches sympy's ``ln`` name used by circuit models
ln = log

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
    
