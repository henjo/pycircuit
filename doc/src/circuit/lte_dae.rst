=================================================
Local Truncation Error (LTE) for DAE Systems
=================================================

Overview
========

Adaptive time-stepping relies on calculating the Local Truncation Error (LTE) to dynamically adjust the step size :math:`h`. However, traditional textbook LTE formulas were derived for ordinary differential equations (ODEs) of the form :math:`dx/dt = f(x)`. 

Circuit simulators, including PyCircuit, actually solve Differential Algebraic Equations (DAEs) of the form:

.. math::
    \frac{d}{dt} q(x(t)) + f(x(t)) + b(t) = 0

where :math:`q(x(t))` represents nonlinear charge or flux. Applying standard ODE LTE formulas to DAEs is merely an approximation that can lead to accuracy and efficiency losses in stiff, highly non-linear circuits.

DAE-Specific LTE Formulation
============================

Following the derivation in *"An Efficient Time Step Control Method in Transient Simulation for DAE System"* (Yao et al.), PyCircuit implements the rigorous LTE formulas specific to DAE systems. 

By applying Taylor expansion with a Lagrange remainder directly to the charge vector :math:`q(x(t))` and substituting it back into the Linear Multi-Step (LMS) equations, we obtain exact analytical formulas for DAE truncation errors.

For example, the DAE LTE for the **Trapezoidal Method** is given by:

.. math::
    \epsilon_T = -\frac{1}{12} (\dot{q}_{x_n} + 0.5 h \ddot{f}_{x_n})^{-1} (g(x_n, t_n) - 2g(x_{n-1}, t_{n-1}) + g(x_{n-2}, t_{n-2})) h

This approach fundamentally ensures that PyCircuit calculates the true truncation error of the charge/flux integration, rather than an ODE approximation.

Impact on Simulation
====================

Using the DAE-specific LTE formulation provides two major benefits:
1. **Accuracy**: Time-step control mathematically matches the true physical dynamics of the non-linear charge nodes.
2. **Efficiency**: It prevents the solver from taking excessively small steps (which happens when ODE approximations break down during fast nonlinear transitions), ultimately speeding up the simulation without sacrificing precision.
