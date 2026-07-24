import time
import warnings
# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

from pycircuit.circuit.circuit import SubCircuit, Node
from pycircuit.circuit.elements import R, VS, Diode
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.toolkit import numeric
import pycircuit.circuit._jaxtoolkit as jax_tk

def build_massive_diode_bridge(num_branches=5000):
    """
    Builds a circuit with a voltage source and `num_branches` parallel diode branches.
    Each branch has a Resistor and a Diode in series.
    This creates a massive non-linear circuit to stress the Newton-Raphson evaluator.
    """
    c = SubCircuit()
    n_in = c.add_node('in')
    gnd = Node('gnd', isglobal=True)
    
    # Input voltage
    c['V1'] = VS(n_in, gnd, v=5.0)
    
    # Massive parallel diode bank
    for i in range(num_branches):
        n_mid = c.add_node(f'mid_{i}')
        # Series Resistor
        c[f'R_{i}'] = R(n_in, n_mid, r=100.0)
        # Non-linear Diode
        c[f'D_{i}'] = Diode(n_mid, gnd, IS=1e-14)
        
    return c

if __name__ == "__main__":
    NUM_BRANCHES = 5000
    print(f"Building a massive circuit with {NUM_BRANCHES} diodes and {NUM_BRANCHES} resistors...")
    circuit = build_massive_diode_bridge(NUM_BRANCHES)
    
    # ---------------------------------------------------------
    # 1. LEGACY PYTHON FOR-LOOP EVALUATION
    # ---------------------------------------------------------
    print("\n--- Running Legacy Python Fallback (numeric toolkit) ---")
    dc_legacy = DC(circuit, refnode=Node('gnd', isglobal=True), toolkit=numeric)
    
    start_time = time.time()
    # Solve DC operating point
    res_legacy = dc_legacy.solve()
    legacy_duration = time.time() - start_time
    print(f"Legacy evaluation completed in: {legacy_duration:.3f} seconds.")
    
    # ---------------------------------------------------------
    # 2. JAX VECTORIZED BATCHED EVALUATION
    # ---------------------------------------------------------
    print("\n--- Running JAX Vectorized Engine (_jaxtoolkit) ---")
    # Tell the new JAX toolkit to warm up and JIT compile the batched functions
    dc_jax = DC(circuit, refnode=Node('gnd', isglobal=True), toolkit=jax_tk)
    
    # The first run will include JAX JIT compilation overhead
    print("JIT Compiling and Solving (First run)...")
    start_time = time.time()
    res_jax = dc_jax.solve()
    jax_first_run = time.time() - start_time
    print(f"JAX (JIT + Solve) completed in: {jax_first_run:.3f} seconds.")
    
    # Second run is purely executing the compiled GPU/CPU kernel
    print("Solving from cached compiled kernel (Second run)...")
    start_time = time.time()
    res_jax2 = dc_jax.solve()
    jax_cached_run = time.time() - start_time
    print(f"JAX (Cached Solve) completed in: {jax_cached_run:.3f} seconds.")
    
    # ---------------------------------------------------------
    # 3. RESULTS & COMPARISON
    # ---------------------------------------------------------
    print("\n=== PERFORMANCE RESULTS ===")
    speedup = legacy_duration / jax_cached_run
    print(f"Speedup vs Legacy: {speedup:.1f}x Faster!")
    
    # Verify exact numerical equivalence
    diff = sum(abs(res_legacy.x - res_jax2.x))
    print(f"Numerical difference between engines: {diff:.2e}")
    if diff < 1e-6:
        print("✅ SUCCESS: Both engines computed the exact same mathematical result!")
    else:
        print("❌ WARNING: The engines produced different results.")
