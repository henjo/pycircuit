import pytest
import numpy as np
from pycircuit.circuit.toolkit import JAXToolkit, NumericToolkit
from pycircuit.circuit import _jaxtoolkit, _numeric
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import Diode, IS

# Instantiate toolkits
jax_tk = JAXToolkit(_jaxtoolkit)
num_tk = NumericToolkit(_numeric)

class JAXDiode(Diode):
    """A diode that uses JAX autodiff instead of manual G(x) and i(x)"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.toolkit = jax_tk

    def eval_i(self, x, epar):
        # We only output one branch (current branch)
        VD = x[0] - x[1] # plus node - minus node
        VT = self.toolkit.kboltzmann * epar.get('T') / self.toolkit.qelectron
        IS = self.ipar.get('IS')
        
        I = IS * (self.toolkit.exp(VD/VT) - 1)
        
        return self.toolkit.array([I, -I])

def test_jax_autodiff_diode():
    # Regular Diode
    d = Diode(0, 1)
    d.toolkit = num_tk
    
    # JAX Diode
    d_jax = JAXDiode(0, 1)
    # generate the eval_i_and_G
    jax_tk.generate_eval_i_and_G(d_jax)
    
    # Test vector
    x = np.array([0.7, 0.0]) # 0.7V forward bias
    
    # Manual evaluation
    manual_i = d.i(x, defaultepar)
    manual_G = d.G(x, defaultepar)
    
    # JAX evaluation
    jax_x = jax_tk.array(x)
    jax_i, jax_G = d_jax.eval_i_and_G(jax_x, defaultepar)
    
    # Check that outputs match
    np.testing.assert_array_almost_equal(manual_i, np.array(jax_i))
    np.testing.assert_array_almost_equal(manual_G, np.array(jax_G))
