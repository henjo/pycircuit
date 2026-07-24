"""
JAX Transient Monte Carlo Analysis
==================================
This example demonstrates a Monte Carlo analysis of a resistive voltage divider.
We introduce a ±5% standard deviation to the resistor values and run 
1,000 simulations to observe the statistical distribution of the output voltage.
Due to JAX's JIT compilation, the solver graph is cached on the first run,
making subsequent simulations run in microseconds.
"""
from pycircuit.circuit.elements import SubCircuit
from pycircuit.circuit import R, VS, gnd
from pycircuit.circuit.jaxtransient import JAXTransient
import pycircuit.circuit._jaxtoolkit as jax_toolkit
import matplotlib.pyplot as plt
import numpy as np
import time

def run_monte_carlo():
    from pycircuit.circuit import circuit
    circuit.default_toolkit = jax_toolkit
    
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    
    # DC input of 1V
    cir['V1'] = VS('in', gnd, v=1.0)
    cir['R1'] = R('in', 'out', r=1.0) # Nominal 1 Ohm
    cir['R2'] = R('out', gnd, r=1.0)  # Nominal 1 Ohm
    
    N_sims = 1000
    std_dev = 0.05 # 5% standard deviation
    
    np.random.seed(42)
    # Generate normally distributed resistor values
    r1_vals = np.random.normal(1.0, std_dev, N_sims)
    r2_vals = np.random.normal(1.0, std_dev, N_sims)
    
    output_voltages = []
    
    print(f"Running {N_sims} Monte Carlo simulations...")
    start_time = time.time()
    
    for r1, r2 in zip(r1_vals, r2_vals):
        # Update component parameters
        cir['R1'].iparv.r = r1
        cir['R2'].iparv.r = r2
        
        # Since this is a DC circuit, solving a very short transient is sufficient
        res = JAXTransient(cir).solve(gnd, tend=1e-5, timestep=1e-5)
        
        # Extract the steady-state output voltage
        v_out = res.x[cir.get_node_index('out')][-1]
        output_voltages.append(v_out)
        
    execution_time = time.time() - start_time
    print(f"Completed in {execution_time:.2f} seconds ({execution_time/N_sims*1000:.2f} ms/sim).")
    
    # Plotting Histogram
    plt.figure(figsize=(8, 6))
    plt.hist(output_voltages, bins=30, color='skyblue', edgecolor='black')
    
    mean_v = np.mean(output_voltages)
    std_v = np.std(output_voltages)
    
    plt.axvline(mean_v, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_v:.3f}V')
    plt.axvline(mean_v + std_v, color='green', linestyle='dashed', linewidth=2, label=f'+1 Std: {mean_v+std_v:.3f}V')
    plt.axvline(mean_v - std_v, color='green', linestyle='dashed', linewidth=2, label=f'-1 Std: {mean_v-std_v:.3f}V')
    
    plt.title('Monte Carlo Analysis: Voltage Divider Output (N=1000)')
    plt.xlabel('Output Voltage (V)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True, alpha=0.5)
    plt.savefig('jax_monte_carlo.png')

if __name__ == '__main__':
    run_monte_carlo()
