=================================================
Local Truncation Error (LTE) for DAE Systems
=================================================

Overview
========

Adaptive time-stepping relies on estimating the Local Truncation Error (LTE) to
dynamically adjust the step size :math:`h`.  Traditional textbook LTE formulas
were derived for ordinary differential equations (ODEs) of the form
:math:`dx/dt = f(x)`.

Circuit simulators, including PyCircuit, actually solve Differential Algebraic
Equations (DAEs) of the form

.. math::
    \frac{d}{dt} q(x(t)) + f(x(t)) + b(t) = 0

where :math:`q(x(t))` represents nonlinear charge or flux.  Applying standard
ODE LTE formulas to DAEs is only an approximation, and can lose accuracy and
efficiency in stiff, highly non-linear circuits.

Two independent layers
======================

It is worth separating two things that are easy to conflate:

1. **The LTE estimate** -- *how* the truncation error :math:`\varepsilon` is
   computed.  This is what the DAE formula below changes.
2. **The step-size selector** -- given :math:`\varepsilon` and a target
   tolerance :math:`\tau`, *how* the next step is chosen, e.g. the elementary
   :math:`h_{new} = (\tau/\varepsilon)^{1/(p+1)}\,h`.

These are orthogonal.  Selecting the DAE LTE formula changes layer 1 only; the
step selector is unchanged.

DAE-specific LTE formulation
============================

Following *"An Efficient Time Step Control Method in Transient Simulation for
DAE System"* (Yao, Wang and Roychowdhury, ICECS 2014), PyCircuit can use the
rigorous LTE formulas specific to DAE systems.  Expanding the charge vector
:math:`q(x(t))` with a Taylor/Lagrange remainder and substituting into the
Linear Multi-Step (LMS) equations yields closed-form DAE truncation errors, e.g.
for the **Trapezoidal** method

.. math::
    \varepsilon_T = -\tfrac{1}{12}\,(\dot{q}_{x_n} + 0.5\,h\,f_{x_n})^{-1}\,
        \big(g(x_n,t_n) - 2 g(x_{n-1},t_{n-1}) + g(x_{n-2},t_{n-2})\big)\,h ,

where :math:`g \equiv \dot{q} = dq/dt`.

The PyCircuit step controller already forms the solution-space error by applying
the Jacobian inverse :math:`J^{-1}` to the ``compute_lte`` residual, and the
companion Jacobian is :math:`J = f_{x_n} + \alpha_0 \dot{q}_{x_n}`.  The paper's
:math:`(\dot{q}_{x_n} + c\,h\,f_{x_n})^{-1}` factor therefore folds directly into
that :math:`J^{-1}`, and each residual reduces to a difference of the charge
derivative :math:`g`:

.. math::
    \text{Euler/GEAR1:}\quad & \varepsilon = -\tfrac{1}{2}\,(g_n - g_{n-1}) \\
    \text{TRAP:}\quad & \varepsilon = -\tfrac{1}{6}\,(g_n - 2 g_{n-1} + g_{n-2}) \\
    \text{GEAR2:}\quad & \varepsilon = -\tfrac{1}{8}\,\frac{h_1+h_2}{h_1 h_2}\,
        \big(h_2\,g_n - (h_1+h_2)\,g_{n-1} + h_1\,g_{n-2}\big)

with :math:`h_1 = t_n - t_{n-1}` and :math:`h_2 = t_{n-1} - t_{n-2}`.

Selecting the formula
=====================

The formula is chosen with the ``lte_formula`` argument of the
:class:`~pycircuit.circuit.integrator.Integrator` strategy passed to
:class:`~pycircuit.circuit.transient.Transient` as ``integrator=``:

* ``lte_formula='ywr'`` (default) -- the Yao--Wang--Roychowdhury DAE formulas
  above;
* ``lte_formula='classic'`` -- the historical divided-difference estimates.

It applies to both the standard adaptive path and the coupled solver
(``coupled_lte=True``), which share the same LTE estimate.

.. code-block:: python

    from pycircuit.circuit.integrator import Gear2Integrator

    tran = Transient(cir, integrator=Gear2Integrator(lte_formula='ywr'))

Classic vs. YWR: a live comparison
==================================

The table below is **generated when this page is built** by charging an RC
circuit (:math:`\tau = 10\,\mu s`) from zero with each integration method under
both LTE formulas, and comparing the final node voltage against the analytic
solution :math:`v(t) = V\,(1 - e^{-t/\tau})`.

.. exec-rst::

    import numpy as np
    from pycircuit.circuit.elements import VS, R, C
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (EulerIntegrator,
                                              TrapezoidalIntegrator,
                                              Gear2Integrator)

    TAU = 10e-6

    def run(integrator_cls, lte):
        c = SubCircuit()
        c['VS'] = VS(1, gnd, v=10)
        c['R1'] = R(1, 2, r=10)
        c['C1'] = C(2, gnd, c=1e-6)
        tran = Transient(c, integrator=integrator_cls(lte), uic=True)
        res = tran.solve(tend=50e-6, timestep=5e-6, coupled_lte=False)
        t = np.asarray(res.sweep_values, dtype=float)
        v_analytic = 10.0 * (1.0 - np.exp(-t[-1] / TAU))
        return len(t), abs(res.v(2, gnd)[-1] - v_analytic)

    rows = []
    for name, integrator_cls in (('euler', EulerIntegrator),
                                 ('trap', TrapezoidalIntegrator),
                                 ('gear2', Gear2Integrator)):
        nc, ec = run(integrator_cls, 'classic')
        ny, ey = run(integrator_cls, 'ywr')
        rows.append((name, nc, ec, ny, ey))

    print(".. list-table:: RC charging (uic, tau = 10 us, max step 5 us):"
          " steps and absolute error vs analytic")
    print("   :header-rows: 1")
    print("   :widths: 14 12 16 12 16")
    print("")
    print("   * - Method")
    print("     - classic steps")
    print("     - classic error")
    print("     - ywr steps")
    print("     - ywr error")
    for name, nc, ec, ny, ey in rows:
        print("   * - %s" % name)
        print("     - %d" % nc)
        print("     - %.2e" % ec)
        print("     - %d" % ny)
        print("     - %.2e" % ey)

Discussion
==========

The comparison isolates exactly where the DAE formula matters:

* **Backward Euler** is *identical* under both formulas -- the derivation above
  shows the YWR GEAR1 residual reduces to the same
  :math:`-\tfrac{1}{2}(g_n - g_{n-1})` the classic path already used.
* **Trapezoidal** differs only in the variable-step bookkeeping; the two agree
  at uniform steps and stay very close in practice.
* **Gear2** is the substantive case.  The classic estimate used a second
  divided difference of the *charge* :math:`q` (effectively a :math:`q''`
  term), which **under-estimates** the true truncation error.  The controller
  is then fooled into taking too few steps, so accuracy suffers.  The YWR form
  uses the second difference of :math:`g = \dot q` -- the correct third-order
  DAE quantity -- and controls the real error.

On the RC transient above this shows up as Gear2 under ``ywr`` taking
appropriately more steps and reaching an error roughly an order of magnitude
smaller than under ``classic`` -- turning Gear2 from the least accurate of the
three methods into the most accurate, as a properly-controlled second-order
method should be.

Conclusion
==========

Using the DAE-specific LTE formulation provides two benefits:

1. **Accuracy** -- time-step control matches the true physical dynamics of the
   non-linear charge nodes rather than an ODE approximation.
2. **Efficiency** -- it neither over-resolves nor (as the classic Gear2
   estimate did) under-resolves fast non-linear transitions, so the solver
   spends steps where they are actually needed.

Because the correction folds into the existing :math:`J^{-1}` step and only
changes the ``compute_lte`` residual, it is a drop-in option shared by both the
adaptive and coupled solvers, with the step-size selector left untouched.
