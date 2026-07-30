import jax.numpy as jnp
import jax

# Simulate example10.py manually
import pycircuit.circuit._jaxtoolkit as jax_toolkit
from pycircuit.circuit.jaxtransient import JAXTransient
from pycircuit.circuit.elements import SubCircuit
from pycircuit.circuit import Circuit, R, VS, gnd, circuit
circuit.default_toolkit = jax_toolkit
import numpy as np

np.random.seed(42)
N_samples = 100
R_nominal = 1e3
R_std = R_nominal * 0.05
R1_samples = np.random.normal(R_nominal, R_std, N_samples)
R2_samples = np.random.normal(R_nominal, R_std, N_samples)

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')

cir['V1'] = VS('in', gnd, v=1.0)
cir['R1'] = R('in', 'out', r=1e3)
cir['R2'] = R('out', gnd, r=1e3)

override_params = {
    'R': {
        'r': jnp.array([R1_samples, R2_samples]).T
    }
}

solver = JAXTransient(cir)
results = solver.solve_batched(gnd, override_params_tree=override_params, tend=1e-12, timestep=1e-12, uic=True)

v_outs = []
for i in range(10):
    v_out = results[i].v('out')[-1]
    v_outs.append(v_out)
    
print("V_outs:", v_outs)
print("R1 samples:", R1_samples[:10])
print("R2 samples:", R2_samples[:10])

