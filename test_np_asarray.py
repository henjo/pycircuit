import jax
import jax.numpy as jnp
import numpy as np

def f(x):
    return np.asarray(x).reshape(-1)

try:
    jax.jit(f)(jnp.array([1, 2]))
    print("Success")
except Exception as e:
    print("Error:", e)
