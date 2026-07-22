Analyses
========

DC
~~

The DC analysis uses Newton-Raphson iterations to solve for the static operating point. 
It leverages the new :doc:`modular_engine` architecture, using continuation methods like 
Gmin-stepping and Source-stepping natively via the ``NonLinearSolver`` interface.

AC (symbolic and numeric)
~~~~~~~~~~~~~~~~~~~~~~~~~

Transient
~~~~~~~~~

Transient analysis simulates the circuit in the time-domain using the :doc:`modular_engine` 
architecture. It supports multiple numerical integration schemes (Euler, Trapezoidal, Gear2), 
adaptive time-stepping (PI and Integral controllers), and unified non-linear solvers (Standard 
and Schur-coupled) to maximize stability and performance.

Noise (symbolic and numeric)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. table:: Parameter summary

  .. parametertable:: pycircuit.circuit.Noise.parameters


2-port (symbolic and numeric)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Loop gain (symbolic and numeric)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

