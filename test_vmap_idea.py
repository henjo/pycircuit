import jax
import jax.numpy as jnp
from typing import NamedTuple

def circuit_i(x, params):
    # params is a dict {'R': R_array, 'C': C_array}
    return x * params['R']

def solve(params):
    def run_chunk(x):
        return circuit_i(x, params)
    
    return run_chunk(jnp.ones(10))

# params = {'R': jnp.ones(10)}
# batched_params = {'R': jnp.ones((5, 10))}
# solve_batched = jax.vmap(solve, in_axes=(0,))
