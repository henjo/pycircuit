"""
JAX Transient Simple Example
============================
This example demonstrates how to use the fast JAX-accelerated Transient solver
for a simple diode clipping circuit.
"""
from pycircuit.circuit import Circuit, R, VS, Diode, gnd
from pycircuit.circuit.jaxtransient import JAXTransient
from pycircuit.circuit.func import Sin
import pycircuit.circuit._jaxtoolkit as jax_toolkit
import matplotlib.pyplot as plt
import numpy as np

def run_example():
    from pycircuit.circuit import circuit
    circuit.default_toolkit = jax_toolkit
    
    from pycircuit.circuit.elements import SubCircuit
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    
    # 5V 1kHz sine wave
    cir['V1'] = VS('in', gnd, v=0)
    cir['V1'].function = Sin(offset=0, amplitude=5.0, freq=1e3, toolkit=jax_toolkit)
    
    cir['R1'] = R('in', 'out', r=1e3)
    cir['D1'] = Diode('out', gnd)
    
    # Solve 3 periods
    solver = JAXTransient(cir)
    res = solver.solve(gnd, tend=0.003, timestep=1e-5)
    
    # Plot results
    t = np.array(res.sweep_values)
    v_in = np.array(res.x[cir.get_node_index('in')])
    v_out = np.array(res.x[cir.get_node_index('out')])
    
    plt.figure()
    plt.plot(t * 1e3, v_in, label='V_in')
    plt.plot(t * 1e3, v_out, label='V_out')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (V)')
    plt.title('JAX Transient Diode Clipper')
    plt.legend()
    plt.grid(True)
    plt.savefig('jax_simple.png')
    
if __name__ == '__main__':
    run_example()
