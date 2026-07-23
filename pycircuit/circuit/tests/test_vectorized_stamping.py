import pytest
import numpy as np
from pycircuit.circuit.circuit import SubCircuit, Node, gnd
from pycircuit.circuit.elements import R, C, Diode
import scipy.sparse
from pycircuit.circuit.toolkit import numeric as num_tk, sparse_numeric as sparse_tk, jaxtoolkit as jax_tk

def test_dense_vectorized_stamping():
    """
    Test that the vectorized coordinate stamping logic creates the correct dense matrix.
    We build a simple resistor divider circuit and manually check the G matrix indices.
    """
    c = SubCircuit()
    c.toolkit = num_tk
    
    n1 = c.add_node('1')
    n2 = c.add_node('2')
    
    # Add two resistors
    # R1: node 1 to node 2
    c['R1'] = R(n1, n2, r=10.0)
    # R2: node 2 to ground
    c['R2'] = R(n2, gnd, r=20.0)
    
    # 3 nodes: 1, 2, gnd
    x = np.zeros(3)
    G = c.G(x)
    
    # Node indices:
    # 0 = node 1
    # 1 = node 2
    # 2 = gnd
    
    # Check G matrix values.
    # R1 adds 1/10 to (0,0), (1,1), -1/10 to (0,1), (1,0)
    # R2 adds 1/20 to (1,1), (2,2), -1/20 to (1,2), (2,1)
    
    assert np.isclose(G[0, 0], 1/10)
    assert np.isclose(G[0, 1], -1/10)
    assert np.isclose(G[1, 0], -1/10)
    assert np.isclose(G[1, 1], 1/10 + 1/20)
    assert np.isclose(G[1, 2], -1/20)
    assert np.isclose(G[2, 1], -1/20)
    assert np.isclose(G[2, 2], 1/20)

def test_sparse_vectorized_stamping():
    """
    Test that the SparseNumericToolkit triggers `build_sparse` and returns a pure scipy.sparse matrix.
    """
    c = SubCircuit()
    c.toolkit = sparse_tk
    
    n1 = c.add_node('1')
    n2 = c.add_node('2')
    
    c['R1'] = R(n1, n2, r=10.0)
    
    x = np.zeros(3)
    G = c.G(x)
    
    # Should be a sparse matrix
    assert scipy.sparse.issparse(G)
    
    # Convert to dense to check values
    G_dense = G.toarray()
    
    assert np.isclose(G_dense[0, 0], 1/10)
    assert np.isclose(G_dense[0, 1], -1/10)
    assert np.isclose(G_dense[1, 0], -1/10)
    assert np.isclose(G_dense[1, 1], 1/10)

def test_jax_vectorized_stamping():
    """
    Test that the JAX toolkit uses `.at[].add()` correctly for immutable tensor assembly.
    """
    import jax
    c = SubCircuit()
    c.toolkit = jax_tk
    
    n1 = c.add_node('1')
    n2 = c.add_node('2')
    
    # Use JAX-compatible elements
    c['R1'] = R(n1, n2, r=5.0)
    
    x = np.zeros(3)
    
    # JAX matrix construction
    G = c.G(x)
    
    assert hasattr(G, 'at')  # It's a jax.numpy array
    
    assert np.isclose(G[0, 0], 1/5)
    assert np.isclose(G[1, 1], 1/5)
    assert np.isclose(G[0, 1], -1/5)
    assert np.isclose(G[1, 0], -1/5)
