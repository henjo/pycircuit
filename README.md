# pycircuit

Python circuit design and simulation tools.

## Install

```bash
pip install -e ".[dev]"
```

## Batched simulation on the JAX backend

`JAXTransient.solve_batched` runs every lane of a parameter sweep
(corners, Monte Carlo) concurrently in one vmapped, jit-compiled
transient — each lane with its own swept parameters and its own DC
operating point.  Build the circuit under the JAX toolkit, then hand
`solve_batched` an override tree keyed by element class:

```python
from pycircuit.circuit import circuit as circuit_mod, gnd
from pycircuit.circuit.toolkit import jaxtoolkit
circuit_mod.default_toolkit = jaxtoolkit   # BEFORE building the circuit

import jax.numpy as jnp
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin, Diode
from pycircuit.circuit.jaxtransient import JAXTransient

cir = SubCircuit()
cir['vs'] = VSin('in', gnd, va=5.0, freq=1e3)
cir['D'] = Diode('in', 'out')
cir['R'] = R('out', gnd, r=1e3)
cir['C'] = C('out', gnd, c=1e-7)

r_lanes = jnp.logspace(2.5, 4.5, 64).reshape(-1, 1)   # (lanes, instances)
results = JAXTransient(cir).solve_batched(
    gnd, override_params_tree={'R': {'r': r_lanes}},
    tend=2e-3, timestep=1e-5)               # one Result per lane
```

Runnable version: `examples/batched_corner_sweep.py`.

When it pays, measured (laptop RTX A1000, details in
`doc/batched_sweep_260821.md`): the batched call overtakes a Python loop
over the CPU `Transient` from roughly 16 lanes, reaching 22.5× at 512
lanes — and its wall-clock barely grows with lane count.  For a handful
of lanes, or for features the traced loop refuses (custom solver
strategies, device bypass), the CPU `Transient` remains the right tool;
the two backends share parameter vocabulary and defaults and are held to
behavioral agreement by a cross-backend conformance test harness
(`doc/backend_parity_260821.md`).

## Documentation

Built with Sphinx. To build locally:

```bash
pip install -e ".[docs]"
make -C doc html
# output in doc/build/html/

# or with Docker:
# docker build -f docker/docs.Dockerfile -t pycircuit-docs .
# docker run --rm -v "$PWD/doc/build:/src/doc/build" pycircuit-docs
```

Hosted on GitHub Pages: https://henjo.github.io/pycircuit/

## Tests

```bash
pytest
```
