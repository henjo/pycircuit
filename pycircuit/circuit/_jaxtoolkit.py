# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric operations using JAX for Auto-Differentiation and JIT
"""

from .constants import *

import jax
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp
import numpy as np
from jax.numpy import cos, sin, tan, cosh, sinh, tanh, log, exp, pi, \
     ceil, floor, dot, linspace, eye, concatenate, sqrt, real, imag,\
     ones, diff, all, maximum, size, conj

alltrue = all

# natural logarithm; matches sympy's ``ln`` name used by circuit models
ln = log

symbolic = False
numeric = True

ac_u_dtype = np.cdouble

def linearsolver(*args, **kvargs):
    return jnp.linalg.solve(*args, **kvargs)

def linearsolverError(*args, **kvargs):
    # JAX doesn't have a specific LinAlgError, but it might raise ValueError
    return ValueError

def toMatrix(array): 
    return array.astype('cdouble')

def det(x): 
    return jnp.linalg.det(x)

def simplify(x): return x

def zeros(*args, **kvargs): 
    return jnp.zeros(*args, **kvargs)

def array(*args, **kvargs): 
    return jnp.array(*args, **kvargs)

def inv(*args, **kvargs): 
    return jnp.linalg.inv(*args, **kvargs)

def integer(x):
    return int(x)

def complex(x):
    return complex(x)
