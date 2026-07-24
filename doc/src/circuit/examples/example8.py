"""
JAX Transient Frequency Sweep
=============================
This example demonstrates sweeping a parameter in time-domain.
We characterize the frequency response of an RC filter by running 
`JAXTransient` iteratively. Since JAX compiles the solver graph 
on the first run, subsequent iterations for different frequencies
execute at lightning speed.
"""
from pycircuit.circuit import Circuit, R, C, VS, gnd
from pycircuit.circuit.jaxtransient import JAXTransient
from pycircuit.circuit.func import Sin
import pycircuit.circuit._jaxtoolkit as jax_toolkit
import matplotlib.pyplot as plt
import numpy as np

def run_sweep():
    from pycircuit.circuit import circuit
    circuit.default_toolkit = jax_toolkit
    
    from pycircuit.circuit.elements import SubCircuit
    cir = SubCircuit()
    cir.add_node('in')
    cir.add_node('out')
    
    # Base circuit
    cir['V1'] = VS('in', gnd, v=0)
    cir['R1'] = R('in', 'out', r=1e3)
    cir['C1'] = C('out', gnd, c=1e-6) # fc = 1 / (2*pi*R*C) = 159 Hz
    
    # Logarithmic sweep from 10Hz to 10kHz
    freqs = np.logspace(1, 4, 30) 
    
    amplitudes = []
    phases = []
    
    print("Running JAX Transient sweep...")
    for i, f in enumerate(freqs):
        cir['V1'].function = Sin(offset=0, amplitude=1.0, freq=f, toolkit=jax_toolkit)
        
        # Simulate 10 periods
        tend = 10.0 / f
        tstep = 1.0 / (f * 100) # 100 points per period
        
        res = JAXTransient(cir).solve(gnd, tend=tend, timestep=tstep)
        
        # Calculate amplitude and phase using the final period
        t = np.array(res.sweep_values)
        v_in = np.array(res.x[cir.get_node_index('in')])
        v_out = np.array(res.x[cir.get_node_index('out')])
        
        period = 1.0 / f
        mask = t >= (tend - 1.1 * period)
        
        t_steady = t[mask]
        v_in_steady = v_in[mask]
        v_out_steady = v_out[mask]
        
        amp = (np.max(v_out_steady) - np.min(v_out_steady)) / 2.0
        amplitudes.append(amp)
        
        # Rough phase estimate via zero-crossing
        in_zero = t_steady[np.where(np.diff(np.sign(v_in_steady)))[0][-1]]
        out_zero = t_steady[np.where(np.diff(np.sign(v_out_steady)))[0][-1]]
        
        time_diff = out_zero - in_zero
        phase = -time_diff * f * 360.0
        # Normalize to [-180, 180]
        while phase > 180: phase -= 360
        while phase < -180: phase += 360
        phases.append(phase)
        
    # Plotting Bode
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    
    ax1.semilogx(freqs, 20 * np.log10(amplitudes))
    ax1.set_title('RC Filter Bode Plot (Transient Sweep)')
    ax1.set_ylabel('Magnitude (dB)')
    ax1.grid(True, which="both", ls="-")
    
    ax2.semilogx(freqs, phases)
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Phase (deg)')
    ax2.grid(True, which="both", ls="-")
    
    plt.tight_layout()
    plt.savefig('jax_sweep_bode.png')

if __name__ == '__main__':
    run_sweep()
