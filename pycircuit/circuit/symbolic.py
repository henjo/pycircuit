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
import numpy as np
from pycircuit.utilities.param import Parameter

def linearsolver(A, b):
    """Numpy compatible wrapper around sympy.solve_linear_system"""

    A = sympy.Matrix(A.tolist())
    b = sympy.Matrix(b.tolist())

#    A.simplify(); b.simplify()

    res = np.array(A.inverse_ADJ() * b)

    return res.reshape((np.size(res,0),) )

def toMatrix(a):
    return sympy.Matrix(a.tolist())

def det(x):
    return sympy.Matrix(x).det()

def setup_analysis(epar):
    """Code that is run by analyses using this toolkit"""
    epar.append(Parameter('kT', default=sympy.Symbol('kT', real=True, positive=True)))
