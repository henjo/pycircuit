"""Batched corner sweep on the JAX backend.

Runs a 64-lane resistor sweep of a half-wave rectifier (VSin -> Diode ->
R || C) in ONE `solve_batched` call: every lane is simulated concurrently
in a single vmapped, jit-compiled loop, and every lane starts from its
OWN DC operating point (solved per lane, with gmin/gshunt/pseudo-
transient continuation as fallbacks).

When this pays: measured on a laptop RTX A1000, the batched call beats a
Python loop over the CPU `Transient` from ~16 lanes up, reaching 22.5x
at 512 lanes -- and the batched wall barely grows with lane count (see
doc/batched_sweep_260821.md).  For a handful of lanes, the plain CPU
loop is the right tool.

Run:  python examples/batched_corner_sweep.py
"""
import numpy as np

## THE ONE SETUP RULE of the JAX backend: the circuit must be BUILT under
## the JAX toolkit, so every element stamps traceable arrays.  Set the
## default toolkit before constructing elements.
from pycircuit.circuit import circuit as circuit_mod, gnd
from pycircuit.circuit.toolkit import jaxtoolkit
circuit_mod.default_toolkit = jaxtoolkit

import jax.numpy as jnp
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin, Diode
from pycircuit.circuit.jaxtransient import JAXTransient

## The circuit -- built ONCE.  Per-lane parameter values come later,
## through the override tree; the values given here are the defaults for
## anything not swept.
cir = SubCircuit()
cir['vs'] = VSin('in', gnd, va=5.0, freq=1e3)
cir['D'] = Diode('in', 'out')
cir['R'] = R('out', gnd, r=1e3)
cir['C'] = C('out', gnd, c=1e-7)

## The sweep: 64 lanes of R, log-spaced over two decades.  The override
## tree is keyed by element CLASS name, then parameter name; the value
## has shape (lanes, instances-of-that-class) -- one R instance here, so
## (64, 1).  Only classes with pure evaluators are sweepable; an
## unknown/unsweepable key raises instead of being silently ignored.
r_lanes = np.logspace(np.log10(3e2), np.log10(3e4), 64)
tree = {'R': {'r': jnp.asarray(r_lanes).reshape(-1, 1)}}

tran = JAXTransient(cir, reltol=1e-4)
results = tran.solve_batched(gnd, override_params_tree=tree,
                             tend=2e-3, timestep=1e-5)

## `results` is a list of per-lane Result objects, each on its own
## adaptive time grid.  Summarize: the rectified DC level (mean of the
## last period) and the ripple, per corner.
print('%10s %12s %12s' % ('R [ohm]', 'mean [V]', 'ripple [V]'))
for r, res in zip(r_lanes, results):
    t = np.asarray(res.sweep_values, float).ravel()
    v = np.asarray(res.v('out'), float).ravel()
    last = t >= 1e-3                      # the second period
    print('%10.1f %12.4f %12.4f'
          % (r, v[last].mean(), v[last].max() - v[last].min()))
