from pycircuit.circuit import Circuit, R, VS, gnd, circuit
import pycircuit.circuit._jaxtoolkit as jax_toolkit
circuit.default_toolkit = jax_toolkit
from pycircuit.circuit.jaxtransient import JAXTransient
from pycircuit.circuit.elements import SubCircuit
import jax.numpy as jnp
import numpy as np

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')

cir['V1'] = VS('in', gnd, v=1.0)
cir['R1'] = R('in', 'out', r=100.0)
cir['R2'] = R('out', gnd, r=100.0)

override_params = {
    'R': {
        'r': jnp.array([[100.0, 100.0]])
    }
}
solver = JAXTransient(cir)
results = solver.solve_batched(gnd, override_params_tree=override_params, tend=1e-12, timestep=1e-12, uic=True)

print(results[0].x)
