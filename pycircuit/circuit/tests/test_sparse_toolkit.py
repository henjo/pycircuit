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


## ---------------------------------------------------------------------------
## STAGE 7d.  These three are xfail(strict), and that is the finding.
##
## They used to build `SubCircuit()` with NO toolkit, so every circuit was
## assembled by the DEFAULT (numeric) backend and only the analysis object was
## handed `sparse_numeric`.  The matrices reaching the solver were dense numpy
## either way, so the tests passed while comparing `numeric` against `numeric`
## with a different solve call bolted on -- they never exercised the path they
## name.
##
## Built with the toolkit under test, THE SPARSE TOOLKIT DOES NOT WORK AT ALL:
## DC, Transient and AC each fail.  `SparseNumericToolkit` makes `cir.G(x)` a
## `coo_matrix`, and `analysis.remove_row_col` calls `numpy.delete` on it, which
## sees a 0-d object and raises `AxisError`.  Every analysis removes the
## reference node, so nothing can complete.  Verified pre-existing: it fails
## identically at `def248c~1`, before stage 7a touched `remove_row_col`.
##
## `strict=True` deliberately: if the sparse path is ever repaired these turn
## red, which is the opposite of how this defect survived the first time.
## ---------------------------------------------------------------------------
_SPARSE_BROKEN = pytest.mark.xfail(
    strict=True,
    reason='sparse_numeric is non-functional: cir.G(x) is a coo_matrix and '
           'remove_row_col calls numpy.delete on it (AxisError). Pre-existing, '
           'see stage 7d. Use linearsolver=AutoSolver() for a sparse solve.')

def _divider(toolkit):
    """Built WITH the toolkit under test -- see the note in test_sparse_dc_analysis."""
    c = SubCircuit(toolkit=toolkit)
    c['IS1'] = IS(gnd, 1, i=1.0)
    c['R1'] = R(1, 2, r=10.0)
    c['R2'] = R(2, gnd, r=10.0)
    return c


@_SPARSE_BROKEN
def test_sparse_dc_analysis():
    """Verify DC analysis with sparse toolkit matches dense.

    STAGE 7d.  These tests used to build `SubCircuit()` with NO toolkit, so every
    circuit was assembled by the DEFAULT (numeric) backend and only the analysis
    was handed `sparse_numeric`.  The matrices reaching the solver were dense
    numpy either way, so the test passed while never exercising the sparse path it
    names -- it compared `numeric` against `numeric` with a different solve call
    bolted on.  The circuits are now constructed with the toolkit under test.
    """
    res_num = DC(_divider(numeric), toolkit=numeric).solve()
    res_sparse = DC(_divider(sparse_numeric), toolkit=sparse_numeric).solve()
    assert_array_almost_equal(res_num.x, res_sparse.x)

def _rc(toolkit):
    c = SubCircuit(toolkit=toolkit)
    c['IS1'] = IS(gnd, 1, i=1.0)
    c['C1'] = C(1, gnd, c=1e-3)
    c['R1'] = R(1, gnd, r=1.0)
    return c


@_SPARSE_BROKEN
def test_sparse_transient_analysis():
    """Verify Transient analysis with sparse toolkit matches dense.

    Circuits built with the toolkit under test -- see test_sparse_dc_analysis.
    """
    res_num = Transient(_rc(numeric), toolkit=numeric).solve(
        tend=1e-3, timestep=1e-4)
    res_sparse = Transient(_rc(sparse_numeric), toolkit=sparse_numeric).solve(
        tend=1e-3, timestep=1e-4)
    assert_array_almost_equal(res_num.x, res_sparse.x)


@_SPARSE_BROKEN
def test_the_sparse_solver_is_actually_reached():
    """The gate the old tests could not state: is the sparse path RUN at all?

    Every assertion above compares two results, and two results agree just as
    happily when both were produced by the dense solver -- which is exactly what
    was happening.  This counts calls instead, so the test fails if the sparse
    backend is bypassed however plausible the numbers look.
    """
    from pycircuit.circuit import _sparse_numeric
    real = _sparse_numeric.linearsolver
    calls = []

    def counting(*a, **k):
        calls.append(1)
        return real(*a, **k)

    _sparse_numeric.linearsolver = counting
    try:
        Transient(_rc(sparse_numeric), toolkit=sparse_numeric).solve(
            tend=1e-3, timestep=1e-4)
    finally:
        _sparse_numeric.linearsolver = real

    assert calls, 'sparse_numeric.linearsolver was never called'
