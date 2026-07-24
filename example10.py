import numpy as np
import jax.numpy as jnp
import matplotlib.pyplot as plt
from pycircuit.circuit import Circuit, R, VS, gnd, circuit
import pycircuit.circuit._jaxtoolkit as jax_toolkit
circuit.default_toolkit = jax_toolkit
from pycircuit.circuit.jaxtransient import JAXTransient
from pycircuit.circuit.elements import SubCircuit
cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')

# DC input of 1V
cir['V1'] = VS('in', gnd, v=1.0)
cir['R1'] = R('in', 'out', r=100.0)
cir['R2'] = R('out', gnd, r=100.0)

# We want to do Monte Carlo on R1 and R2
# +/- 5% standard deviation normal distribution
N_samples = 10000
np.random.seed(42)

# Generate parameters
R1_samples = np.random.normal(100.0, 5.0, N_samples)
R2_samples = np.random.normal(100.0, 5.0, N_samples)

# Setup override_params_tree for solve_batched
# The tree matches the internal JAX grouping.
# For elements like R, the internal group name is 'R'
override_params = {
    'R': {
        'r': jnp.array([R1_samples, R2_samples]).T  # shape (N_samples, 2)
    }
}

solver = JAXTransient(cir)

# We run a very short transient (e.g. 1 step) to find the DC equivalent output
# Since there are no capacitors, the initial state (even if 0) will jump to DC instantly
results = solver.solve_batched(gnd, override_params_tree=override_params, tend=1e-12, timestep=1e-12, uic=True)

Vout_dc = np.array([res.v('out')[-1] for res in results])

print(f"Mean Vout: {np.mean(Vout_dc):.4f} V")
print(f"Std Vout:  {np.std(Vout_dc):.4f} V")

# Plot histogram
plt.hist(Vout_dc, bins=50, edgecolor='black')
plt.title('Monte Carlo Output Voltage (10,000 runs)')
plt.xlabel('Voltage (V)')
plt.ylabel('Frequency')
plt.savefig('example10_monte_carlo.png')
print("Saved example10_monte_carlo.png")
