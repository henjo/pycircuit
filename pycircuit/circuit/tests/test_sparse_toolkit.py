import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, C, IS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.toolkit import numeric, sparse_numeric

def test_sparse_numeric_linearsolver():
    """Verify that sparse_numeric provides the same linear solve result as numeric."""
    # Simple linear system A*x = b
    # A = [ 2.0 -1.0 ]
    #     [-1.0  2.0 ]
    # b = [ 1.0 ]
    #     [ 0.0 ]
    A = np.array([[2.0, -1.0], [-1.0, 2.0]], dtype=float)
    b = np.array([1.0, 0.0], dtype=float)
    
    x_num = numeric.linearsolver(A, b)
    x_sparse = sparse_numeric.linearsolver(A, b)
    
    assert_array_almost_equal(x_num, x_sparse)

def test_sparse_dc_analysis():
    """Verify DC analysis with sparse toolkit matches dense."""
    c = SubCircuit()
    c['IS1'] = IS(gnd, 1, i=1.0)
    c['R1'] = R(1, 2, r=10.0)
    c['R2'] = R(2, gnd, r=10.0)
    
    # Dense solve
    dc_num = DC(c, toolkit=numeric)
    res_num = dc_num.solve()
    
    # Sparse solve
    dc_sparse = DC(c, toolkit=sparse_numeric)
    res_sparse = dc_sparse.solve()
    
    assert_array_almost_equal(res_num.x, res_sparse.x)

def test_sparse_transient_analysis():
    """Verify Transient analysis with sparse toolkit matches dense."""
    c = SubCircuit()
    c['IS1'] = IS(gnd, 1, i=1.0)
    c['C1'] = C(1, gnd, c=1e-3)
    c['R1'] = R(1, gnd, r=1.0)
    
    # Dense solve
    tran_num = Transient(c, toolkit=numeric)
    res_num = tran_num.solve(tend=1e-3, timestep=1e-4)
    
    # Sparse solve
    tran_sparse = Transient(c, toolkit=sparse_numeric)
    res_sparse = tran_sparse.solve(tend=1e-3, timestep=1e-4)
    
    assert_array_almost_equal(res_num.x, res_sparse.x)
