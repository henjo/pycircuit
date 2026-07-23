# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of sparse numeric operations that can be used as a toolkit for
Analysis objects

The module is based on `scipy.sparse <http://scipy.org>`_ and wraps
`pycircuit.circuit._numeric`.

"""

from ._numeric import *
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, MatrixRankWarning
import warnings

def linearsolver(*args, **kvargs):
    A = args[0]
    b = args[1]
    
    A_sparse = csc_matrix(A)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", MatrixRankWarning)
        
        try:
            res = spsolve(A_sparse, b, **kvargs)
        except RuntimeError as e:
            if "singular" in str(e).lower():
                raise np.linalg.LinAlgError("Singular matrix")
            raise e
            
        if len(w) > 0 and issubclass(w[-1].category, MatrixRankWarning):
            raise np.linalg.LinAlgError("Singular matrix")
            
        if np.any(np.isnan(res)):
            raise np.linalg.LinAlgError("Singular matrix")
            
        # spsolve returns a rank-1 array if b is rank-1.
        # But we ensure it's compatible with what numeric returns
        return res

numeric = True
sparse = True
