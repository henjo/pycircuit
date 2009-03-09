Theory
------

Solving
```````

pycircuit.circuit uses the [MNA] (Modified Nodal Analysis) formulation to 
solve the voltages and currents of an electric circuit.

For a circuit with N nodes and M independent voltage sources or inductors,
the unknown quantities in the MNA method can be written as:

.. math::
  x = \left( 
        \begin{array}{c}
        v\\
        j 
        \end{array}
      \right)

where :math:`v` is a vector of length N that contains the node
voltages and :math:`j` is a vector of branch currents through
independent voltage sources and inductors.

The behaviour of a non-linear lumped circuit with energy storage elements 
such as inductors and capacitors can be found from the solution of a system of
differential equations:

.. math::
  \frac{d}{dt} q(v(t)) + i(v(t)) + u(t) = 0
  :label: diffeqnonlin

Where :math:`t` denotes time, :math:`v(t)` is a vector of node voltages and terminal currents. :math:`q(\cdot)` is a function that maps :math:`v(t)` to sums of capacitive charge or inductive fluxes at a node, :math:`i(\cdot)` is a function that maps the :math:`x(t)` vector to a sum of currents at a node and finally 
:math:`u(t)` is the input to the system.

If the circuit is linear or linearized around the operating point, equation :eq:`diffeqnonlin` can be rewritten in matrix form as:

.. math::
  C \frac{d}{dt} x(t) + G x(t) + u(t) = 0 
  :label: eqlinearcircuit

where the matrix elements are given by:

.. math::
  C_{j, k} = \frac{\partial q_j}{\partial x_k} \\
  G_{j, k} = \frac{\partial i_j}{\partial x_k}

AC analysis
```````````

In an AC analysis the stimuli the circuit is a complex sinusoid, that is:

.. math::
  u(t) = A e^{j (\omega t + \phi)}

The solution can be found by applying the Laplace transform to equation
:eq:`eqlinearcircuit`:

.. math::
  \begin{array}{ll}
    \LARGE{\mathcal{L}} \left\{ C \frac{d}{d t} x(t) + G x(t) + u(t)
    \right\} = \large{\LARGE{\mathcal{L}}} \left\{ 0 \right\} & \Rightarrow \\
    s C X(s) + G X(s) + U(s) = 0 & \Rightarrow\\
    (G + s C) X(s) + U(s) = 0 & \Rightarrow\\
    X(s) = -(G + s C)^{- 1} U(s)
  \end{array}


.. [MNA]: http://en.wikipedia.org/wiki/Modified_nodal_analysis
