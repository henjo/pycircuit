# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of symbolic operations that can be used as a toolkit for
Analysis objects.

The module is based on the python CAS `sympy <http://sympy.org>`_.

"""

import sympy
from sympy import cos, sin, tan, exp, pi, simplify, floor
from sympy import oo as inf, ceiling as ceil
from numpy import dot, eye, linspace
import numpy as np
import types
from pycircuit.utilities.param import Parameter

symbolic = True

ac_u_dtype = np.object


def linearsolver(A, b):
    A,subst_dict = dummy_var_matrix(A)
    b = sympy.Matrix(b.tolist())

#    A.simplify(); b.simplify()

    res = np.array((A.inverse_ADJ() * b).subs(subst_dict))

    return res.reshape((np.size(res,0),) )

def toMatrix(a):
    return sympy.Matrix(a.tolist())

def det(A):
    A, subst_dict = dummy_var_matrix(A)
    return A.det().subs(subst_dict)

def cofactor(x, i, j):
    return sympy.Matrix(x).cofactor(i, j)

def setup_analysis(epar):
    """Code that is run by analyses using this toolkit"""
    epar.append(Parameter('kT', default=sympy.Symbol('kT', real=True, positive=True)))

def zeros(shape, dtype=object):
    return np.zeros(shape, dtype=dtype)

def array(x, dtype=object):
    return np.array(x, dtype=dtype)

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
    Aprime = sympy.Matrix(A.lines, A.cols, elem)
    return Aprime, subst_dict


