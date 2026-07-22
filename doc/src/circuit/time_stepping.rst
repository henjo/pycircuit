================================================
Transient Analysis & Time-Stepping Architecture
================================================

This document explains the mathematical architecture behind PyCircuit's transient simulation engine, specifically detailing the integrators, nonlinear solvers, and adaptive time-stepping methods.

Integrators (Differentiators)
=============================

PyCircuit solves Differential Algebraic Equations (DAEs) of the form:

.. math::
    f(v(t)) + \frac{d}{dt} q(v(t)) + u(t) = 0

To solve this numerically, the continuous time derivative must be replaced by a discrete-time approximation. This is the job of the **Integrator** (often referred to as a differentiator in MNA context).

PyCircuit supports two primary implicit integration methods:

* **Trapezoidal Rule**: A 2nd-order method that uses a linear interpolation of the derivative. It is highly accurate and preserves energy in non-dissipative circuits (like LC tanks) but can suffer from numerical "ringing" (oscillations) during sharp discontinuous transitions.
* **Gear-2 (BDF2)**: A 2nd-order Backward Differentiation Formula. It is strictly A-stable and heavily damps numerical ringing, making it the industry standard for stiff circuits with sharp edges, though it slightly dampens physical oscillations in lossless circuits.

Nonlinear Solvers
=================

At every discrete time step :math:`t_n`, the discretized DAE becomes a system of nonlinear algebraic equations:

.. math::
    F(x_n) = 0

To find the circuit state :math:`x_n` (voltages and branch currents), PyCircuit uses the **Newton-Raphson (NR)** method.

* The solver computes the Jacobian matrix :math:`J = \partial F / \partial x`.
* It iteratively updates the state: :math:`\Delta x = -J^{-1} F(x)`.
* The NR loop terminates when both the absolute and relative changes in :math:`x` fall below user-defined tolerances (``vabstol``, ``iabstol``, ``reltol``).

Adaptive Time-Stepping Methods
==============================

Choosing the right step size :math:`h = t_n - t_{n-1}` is critical. Small steps are needed during fast transitions to maintain accuracy, while large steps are needed during inactive periods for simulation speed.

Option B: Standard Predictive Time-Stepping (Legacy)
----------------------------------------------------

This is the traditional "trial and error" approach used by SPICE:

1. The solver predicts a step size :math:`h` based on the Local Truncation Error (LTE) of the *previous* time steps.
2. It solves the nonlinear NR system for the current state :math:`x_n` using :math:`h`.
3. After convergence, it calculates the actual LTE for the *current* step.
4. If the LTE exceeds the error tolerance (``TRTOL``), the step is **rejected** (backed up). :math:`h` is reduced, and the solver tries again.

**Drawbacks**: For highly stiff circuits, the predicted :math:`h` is often overly optimistic, leading to frequent backups and wasted matrix inversions.

Option A: Coupled Schur Complement Time-Stepping (New)
------------------------------------------------------

This method (based on G. Peter Fang's "A New Time-Stepping Method for Circuit Simulation") eliminates rejected time steps entirely by treating the timestep :math:`h` as an independent variable to be solved *simultaneously* with the circuit state :math:`x`.

The solver couples the LTE error equation :math:`E(x_n, h) = 0` directly into the DAE system, forming an augmented :math:`(N+1) \times (N+1)` system:

1. **Approximate Newton Update**: To avoid solving a massive :math:`(N+1)` matrix, PyCircuit exploits the analytical relationship between LTE and :math:`h` (where LTE :math:`\propto h^3` for Trapezoidal). 
2. **Analytical Gradient** (:math:`E_h`): The sensitivity of the error to the timestep is calculated analytically as :math:`E_h = \frac{p \cdot (E + TRTOL)}{h}`, which provides a perfectly smooth, mathematically stable Newton update for :math:`h`.
3. **Golden Window** (:math:`\gamma` bounds): Following Section 3.3 of Fang's paper, the solver defines an acceptable error window (:math:`0.7\tau \le \epsilon \le 3.0\tau`). If the current error falls inside this window, the solver accepts the step without trying to over-optimize :math:`h`, saving thousands of redundant iterations.

This approach guarantees that the solver scales :math:`h` aggressively upward when the circuit is quiet, and shrinks :math:`h` stably without getting trapped by numerical noise during severe stiffness.
