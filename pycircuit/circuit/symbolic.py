# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of symbolic operations that can be used as a toolkit for
Analysis objects.

The module is based on the python CAS `sympy <http://sympy.org>`_.

"""

import numpy as np
import sympy
from sympy import cos, sin, tan, exp, pi, simplify, floor
from sympy import oo as inf, ceiling as ceil
import types
from pycircuit.utilities.param import Parameter
from pycircuit.circuit.constants_sympy import kboltzmann, eps0, epsRSi, epsRSiO2, qelectron

symbolic = True

ac_u_dtype = object

# def linearsolver(A, b):
#     A,subst_dict = dummy_var_matrix(A)
#     b = sympy.Matrix(b.tolist())

# #    A.simplify(); b.simplify()

#     if A.shape == (1,1):
#         return np.array([(b[0] / A[0,0]).subs(subst_dict)])
#     else:
#         res = np.array((A.inverse_ADJ() * b).subs(subst_dict))

#     return res.reshape((np.size(res,0),) )

def linearsolver(A, b):
    A = sympy.Matrix(A)
    b = sympy.Matrix(b.tolist())

#    A.simplify(); b.simplify()

    if A.shape == (1,1):
        return np.array([(b[0] / A[0,0])])
    else:
        #res = np.array((A.inverse_ADJ() * b))
        res = np.array((A.LUsolve(b)))

    return res.reshape((np.size(res,0),) )

def linearsolverError(*args, **kvargs):
    return np.linalg.LinAlgError

def toMatrix(a):
    return sympy.Matrix(a.tolist())

def det(A):
    A, subst_dict = dummy_var_matrix(A)
    return A.det().subs(subst_dict)

def cofactor(x, i, j):
    return sympy.Matrix(x).cofactor(i, j)

def setup_analysis(epar):
    """Code that is run by analyses using this toolkit"""
    epar['T'] = Parameter('T', default=sympy.Symbol('T', real=True, positive=True))

def zeros(shape, dtype=None):
    return np.zeros(shape, dtype=object)

def nonzero(x):
    """Return a list of non-zero indices of array x"""
    return np.where(x != 0)[0]

def dummy_var_matrix(A):
    """Substitute non-zero elements with symbols
    
    Returns the new matrix and a dictionary that can be used
    with the subs method to go back to the original matrix.
    
    This can improve performance considerably of the derminant function in
    Sympy.
    """

    A = sympy.Matrix(A)
    subst_dict = {}
    def elem(i,j):
        if A[i,j] != 0:
            sym = sympy.Symbol('a%d%d'%(i,j))
            subst_dict[sym] = A[i,j]
            return sym
        else:
            return 0
    Aprime = sympy.Matrix(A.rows, A.cols, elem)
    return Aprime, subst_dict

def array(*args,**kvargs):
    return np.array(*args,**kvargs)

def dot(*args,**kvargs):
    return np.dot(*args,**kvargs)

def delete(*args,**kvargs):
    return np.delete(*args,**kvargs)

def eye(*args,**kvargs):
    return np.eye(*args,**kvargs)

def linspace(*args,**kvargs):
    return np.linspace(*args,**kvargs)

def ones(*args,**kvargs):
    return np.ones(*args,**kvargs)

def concatenate(*args,**kvargs):
    return np.concatenate(*args,**kvargs)

def imag(*args,**kvargs):
    return sympy.im(*args,**kvargs)

def conj(xarray):
    conjlist = []
    for x in xarray:
        conjlist.append(sympy.conjugate(x))
    return array(conjlist)

def size(*args,**kvargs):
    return np.size(*args,**kvargs)

def integer(x):
    return sympy.Integer(x)
