# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Module of numeric operations using JAX for Auto-Differentiation and JIT
"""

from .constants import *

import jax
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp
import numpy as np
from jax.numpy import cos, sin, tan, cosh, sinh, tanh, log, exp, pi, \
     ceil, floor, dot, linspace, eye, concatenate, sqrt, real, imag,\
     ones, diff, all, maximum, size, conj, inf

alltrue = all

# natural logarithm; matches sympy's ``ln`` name used by circuit models
ln = log

symbolic = False
numeric = True

ac_u_dtype = np.cdouble

def linearsolver(*args, **kvargs):
    return jnp.linalg.solve(*args, **kvargs)

def linearsolverError(*args, **kvargs):
    # JAX doesn't have a specific LinAlgError, but it might raise ValueError
    return ValueError

def toMatrix(array): 
    return array.astype('cdouble')

def det(x): 
    return jnp.linalg.det(x)

def simplify(x): return x

def zeros(*args, **kvargs): 
    return jnp.zeros(*args, **kvargs)

def array(*args, **kvargs): 
    return jnp.array(*args, **kvargs)

def inv(*args, **kvargs): 
    return jnp.linalg.inv(*args, **kvargs)

def integer(x):
    return int(x)

def complex(x):
    return complex(x)

def sum(*args, **kvargs):
    return jnp.sum(*args, **kvargs)

def reshape(a, newshape):
    return jnp.reshape(a, newshape)

def add_at(a, indices, b):
    # indices is a tuple of (rows, cols)
    return a.at[indices].add(b)

def delete(arr, obj, axis=None):
    if isinstance(obj, list):
        obj = jnp.array(obj)
    return jnp.delete(arr, obj, axis=axis)

def insert(arr, obj, values, axis=None):
    if isinstance(obj, list):
        obj = jnp.array(obj)
    return jnp.insert(arr, obj, values, axis=axis)




where = jnp.where
logical_and = jnp.logical_and
logical_or = jnp.logical_or

jax = True

def generate_batched_eval(element_cls, method='i'):
    import jax
    import pycircuit.circuit._jaxtoolkit as toolkit
    
    cache_attr = f'_jax_batched_{method}'
    if not hasattr(element_cls, cache_attr):
        setattr(element_cls, cache_attr, {})
        
    cache = getattr(element_cls, cache_attr)
    
    def batched_eval_func(X_batch, params_batch, epar):
        if hasattr(epar, '_values'):
            key = frozenset(epar._values.items())
        else:
            key = id(epar)
            
        if key not in cache:
            @jax.jit
            def compiled_func(X_b, p_b):
                def single_eval(x, p):
                    if method == 'i':
                        return element_cls.eval_i_pure(x, p, epar, toolkit)
                    elif method == 'q':
                        if hasattr(element_cls, 'eval_q_pure'):
                            return element_cls.eval_q_pure(x, p, epar, toolkit)
                        else:
                            import jax.numpy as jnp
                            return jnp.zeros_like(element_cls.eval_i_pure(x, p, epar, toolkit))
                
                # --- HIGH-PERFORMANCE GPU/CPU VECTORIZATION ---
                # Instead of evaluating thousands of Diodes/Transistors in a Python `for` loop,
                # we use JAX's `vmap` (Vectorizing Map). It takes our `single_eval` function 
                # (which operates on a single device) and maps it across a massive stacked array 
                # of inputs (X_b) and parameters (p_b). 
                # 
                # Simultaneously, `jax.jacfwd` performs Forward-Mode Automatic Differentiation 
                # to EXACTLY calculate the Jacobian (the Conductance / Capacitance matrix) for 
                # every single device without manually writing calculus formulas.
                # 
                # JAX compiles this entire block via `@jax.jit` into a highly optimized 
                # XLA executable that executes simultaneously in C++/GPU space.
                val_batch = jax.vmap(single_eval)(X_b, p_b)
                jac_batch = jax.vmap(jax.jacfwd(single_eval))(X_b, p_b)
                return val_batch, jac_batch
                
            cache[key] = compiled_func
            
        return cache[key](X_batch, params_batch)
        
    return batched_eval_func
    
