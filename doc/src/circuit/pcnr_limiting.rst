===============================================================
Predictor-Corrector Newton-Raphson (PCNR) Limiting
===============================================================

Overview
========

In non-linear circuit simulation, Newton-Raphson (NR) convergence can easily fail due to exponential device relationships (e.g., P-N junctions in Diodes and BJTs). If a node voltage jumps too drastically in a single iteration, evaluating :math:`e^{V/V_T}` can cause numeric overflows or wildly unstable Jacobians.

Traditionally, SPICE simulators solve this using a practice called "limiting" (e.g., the ``pnjlim`` algorithm). Standard limiting forces devices to maintain an internal, hidden "junction voltage" history between iterations to artificially restrict voltage swings. 

**Problems with Traditional Limiting**:
* The Jacobian matrix becomes mathematically inconsistent because devices evaluate currents at a different voltage than the rest of the MNA matrix expects.
* It breaks modularity by forcing devices to track simulation history and iteration states.

PCNR: A Consistent Replacement
==============================

PyCircuit implements **Predictor-Corrector Newton-Raphson (PCNR)**, based on *"Predictor/Corrector Newton-Raphson (PCNR): A Simple, Flexible, Scalable, Modular, and Consistent Replacement for Limiting in Circuit Simulation"* (Aadithya et al.).

PCNR replaces traditional stateful limiting with a perfectly consistent mathematical abstraction.

1. **Explicit Unknowns**: Instead of hiding limited junction voltages inside the device, PCNR elevates every limited quantity to a formal, explicit unknown variable in the circuit's MNA system. 
2. **Predictor Phase**: The NR solver computes an unconstrained mathematical update for the entire circuit using standard MNA matrices.
3. **Corrector Phase**: Each limited device explicitly "limits" its own assigned variables, and these limits are globally passed to the next iteration.

Schur Complement Acceleration
=============================

Naively adding explicit unknowns for every limited device would vastly increase the size of the MNA matrix, slowing down the sparse matrix solvers (``Ax=b``). 

To prevent this overhead, PyCircuit uses a **Schur Complement** reduction technique. Because the block sub-matrix associated with the new limited variables is strictly the Identity matrix, the extra variables can be analytically eliminated during the matrix solve. 

This guarantees that PyCircuit executes the PCNR limiting math perfectly while maintaining an MNA matrix size identical to that of traditional SPICE, achieving perfect mathematical consistency without any performance penalty.
