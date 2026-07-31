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

    import warnings
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
        ## Counted rather than printed: a force-accept is an accepted step whose
        ## truncation error nothing bounds, so it belongs in the table next to the
        ## error it explains -- not in the build log, where it reads as a defect
        ## in the build.  See "The step-size ratio, and what bounds it" below.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            res = tran.solve(tend=50e-6, timestep=5e-6, coupled_lte=False)
        forced = sum(1 for w in caught
                     if issubclass(w.category, RuntimeWarning)
                     and 'truncation error' in str(w.message))
        t = np.asarray(res.sweep_values, dtype=float)
        v_analytic = 10.0 * (1.0 - np.exp(-t[-1] / TAU))
        return len(t), abs(res.v(2, gnd)[-1] - v_analytic), forced

    rows = []
    for name, integrator_cls in (('euler', EulerIntegrator),
                                 ('trap', TrapezoidalIntegrator),
                                 ('gear2', Gear2Integrator)):
        nc, ec, fc = run(integrator_cls, 'classic')
        ny, ey, fy = run(integrator_cls, 'ywr')
        rows.append((name, nc, ec, fc, ny, ey, fy))

    print(".. list-table:: RC charging (uic, tau = 10 us, max step 5 us):"
          " steps, absolute error vs analytic, and force-accepts")
    print("   :header-rows: 1")
    print("   :widths: 12 11 14 11 11 14 11")
    print("")
    print("   * - Method")
    print("     - classic steps")
    print("     - classic error")
    print("     - classic forced")
    print("     - ywr steps")
    print("     - ywr error")
    print("     - ywr forced")
    for name, nc, ec, fc, ny, ey, fy in rows:
        print("   * - %s" % name)
        print("     - %d" % nc)
        print("     - %.2e" % ec)
        print("     - %d" % fc)
        print("     - %d" % ny)
        print("     - %.2e" % ey)
        print("     - %d" % fy)

The last column counts steps the solver accepted with an unbounded truncation
error, having rejected them three times first.  On this circuit only one
configuration reaches that path, and it is the one whose LTE formula is a
uniform-grid formula being used on a non-uniform grid -- see the discussion of
Trapezoidal below, and the section on the step-size ratio that follows it.

Discussion
==========

The comparison isolates exactly where the DAE formula matters:

* **Backward Euler** is *identical* under both formulas -- the derivation above
  shows the YWR GEAR1 residual reduces to the same
  :math:`-\tfrac{1}{2}(g_n - g_{n-1})` the classic path already used.
* **Trapezoidal** differs only in the variable-step bookkeeping; the two agree
  at uniform steps and stay very close in practice.
* **Gear2** is the substantive case.  Until 2026-07 the classic estimate used a
  second divided difference of the *charge* :math:`q` -- which is a
  :math:`q''` term -- scaled by :math:`h^3`, where BDF-2 requires
  :math:`q'''` scaled by :math:`h^2`.  The YWR form uses the second difference
  of :math:`g = \dot q`, the correct third-order DAE quantity.

.. warning::

   **This page previously understated that defect, and the understatement was
   the more dangerous half.**  It described the consequence as the controller
   being "fooled into taking too few steps, so accuracy suffers", and quantified
   it as "roughly an order of magnitude" of error.  The controller was not
   taking *too few* steps; it was **not controlling the step at all**.

   The mismatched estimate is not a constant factor -- it is wrong by a factor of
   order :math:`h\omega`, which on a 1 MHz signal at nanosecond steps is
   :math:`\sim 10^{-15}`.  With ``err`` that small against an accept threshold of
   1, the growth limiter :math:`\min(2,\,0.9\,(1/\mathrm{err})^{1/3})` saturates
   at 2 on *every* step, so :math:`h` doubles until it reaches ``max_step`` and
   pins there.  Measured consequences: **zero step rejections ever**, and results
   *bit-identical* across a :math:`10^3` change in ``abstol`` and ``reltol``.
   Tolerances would have to be some :math:`3\times10^{10}` times tighter to
   provoke a single rejection.

   A blind controller silently delivers whatever accuracy ``timestep`` happens to
   buy, and tightening tolerances buys nothing.  That is a different failure from
   "less accurate", and describing it as the latter is why the defect survived
   being documented.

   **Fixed 2026-07-30.**  ``'classic'`` now estimates :math:`q'''` as twice the
   second divided difference of :math:`g`, taken from the companion-current
   history -- the same information YWR uses, since Gear-2 keeps only two past
   charges and a third divided difference of :math:`q` is therefore unavailable.
   The estimate is asymptotically exact (ratio to the true truncation error
   :math:`0.9795 \to 0.9988` as :math:`h` halves, against a flat
   :math:`\sim\!10^{-16}` before), and on a stiff reference case the controller
   went from 0 rejections and 0.0% step-count response to ``reltol``, to 4
   rejections and +327.5%.  ``Gear2Integrator`` now also *defaults* to ``'ywr'``.
   See ``doc/transient_repair_plan.md`` for the full gate outcomes.

The step-size ratio, and what bounds it
=======================================

Everything above is about *how large* the truncation error is.  There is a second
constraint that has nothing to do with accuracy: how fast the step size may
**change**.

Variable-step BDF-2 fits a quadratic through :math:`t_n`, :math:`t_{n-1}` and
:math:`t_{n-2}`, and the homogeneous part of that recursion has two roots.  One is
:math:`z = 1`, which is what consistency requires.  The other is *parasitic*, and
its magnitude depends only on the ratio :math:`r = h_n/h_{n-1}`:

.. math::

   z_{\text{parasitic}} = \frac{r^2}{2r + 1}

The method is zero-stable only while that root stays inside the unit disc, i.e.
while :math:`r^2 < 2r+1`, i.e. :math:`r < 1 + \sqrt{2} \approx 2.414214`
(Grigorieff).  Above it, an error already present in the solution is *amplified*
by every subsequent step instead of being forgotten:

.. exec-rst::

    import math

    def parasitic(r):
        return r * r / (2.0 * r + 1.0)

    bound = 1.0 + math.sqrt(2.0)
    rows = [(bound, 'the bound itself'),
            (2.5, 'just past it'),
            (3.0, ''),
            (10.0, 'the escape hatch, before 4b')]

    print(".. list-table:: Parasitic root of variable-step BDF-2, and what it does"
          " over 20 steps")
    print("   :header-rows: 1")
    print("   :widths: 14 18 20 24")
    print("")
    print("   * - ratio r")
    print("     - parasitic root")
    print("     - growth over 20 steps")
    print("     - where this ratio comes from")
    for r, note in rows:
        z = parasitic(r)
        print("   * - %.4f" % r)
        print("     - %.6f" % z)
        print("     - %.3g" % (z ** 20))
        print("     - %s" % (note or "-"))

Two clamps keep a run inside that bound, and they are deliberately different
things:

* ``stepcontroller.MAX_GROWTH_RATIO`` (2.0) is the factor the step controller
  will never exceed when predicting the next step.  It sits *inside* the
  stability bound rather than on it: at exactly :math:`1+\sqrt 2` the parasitic
  root is 1.000000, and a method running permanently at its own stability
  boundary has no margin for the rounding that put it there.
* ``integrator.ZERO_STABILITY_RATIO`` (:math:`1+\sqrt 2`) is a backstop inside
  ``Gear2Integrator.check_order_drop``.  If a ratio above the bound ever reaches
  the integrator, it drops to backward Euler for that step -- order 1 has no
  parasitic root to amplify, so the ratio becomes harmless rather than forbidden.

The backstop is expected to be idle.  Its value is that it makes the bound a
property the integrator enforces for itself, rather than one that happens to hold
because the controller's clamp is smaller.

.. note::

   **The one path that used to bypass both.**  When the truncation error stays
   above tolerance for several successively smaller steps -- which happens near a
   source discontinuity, where the estimate is built on differences of something
   that is not differentiable -- the solver gives up rejecting, accepts the
   converged solution, and moves on.  That escape hatch used to set the next step
   to :math:`10\times` the current one: growing tenfold in answer to an error that
   was already too large, at a parasitic root of 4.76, and silently.

   It now drops the order instead, grows by ``MAX_GROWTH_RATIO`` like every other
   accepted step, and raises a ``RuntimeWarning`` naming the time and the step
   size.  An accepted truncation error that nothing bounds must not be invisible.

The table below is measured when this page is built, on a run chosen because it
*does* reach that escape hatch:

.. exec-rst::

    import warnings
    import numpy as np
    from pycircuit.circuit.elements import R, C, L
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (TrapezoidalIntegrator,
                                              Gear2Integrator,
                                              ZERO_STABILITY_RATIO)
    from pycircuit.circuit.stepcontroller import IntegralController

    class Spy(IntegralController):
        """Records the ACCEPTED step sequence, force-accepted steps included.

        A force-accepted step enters the integrator history like any other, so it
        belongs in the ratio sequence.
        """
        def __init__(self):
            super().__init__()
            self.h = []
            self.forced = 0
            self._consec = 0

        def evaluate_step(self, *args, **kwargs):
            accept, h_next = super().evaluate_step(*args, **kwargs)
            forced = False
            if not accept:
                self._consec += 1
                if self._consec > 3:          # transient.py's MAX_REJECT
                    forced, self._consec = True, 0
                    self.forced += 1
            else:
                self._consec = 0
            if accept or forced:
                self.h.append(kwargs['h_curr'])
            return accept, h_next

    def run(integrator, reltol):
        c = SubCircuit()
        c['C1'] = C(1, gnd, c=1e-6)
        c['R1'] = R(1, 2, r=1.0)
        c['L1'] = L(2, gnd, L=1e-6)
        x0 = np.zeros(c.n)
        x0[c.get_node_index('1')] = 1.0
        x0[c.get_node_index('2')] = 1.0
        tran = Transient(c, integrator=integrator, reltol=reltol)
        spy = Spy()
        tran.step_controller = spy
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            tran.solve(tend=5e-3, timestep=2e-4, x0=x0, coupled_lte=False)
        ratios = [spy.h[i + 1] / spy.h[i] for i in range(len(spy.h) - 1)]
        return len(spy.h), spy.forced, max(ratios), \
            sum(1 for r in ratios if r > ZERO_STABILITY_RATIO)

    ## Both tolerances matter.  Gear2 -- the default -- only reaches the escape
    ## hatch at the LOOSER one, which is where a step is free to grow into
    ## trouble; a sweep of tight tolerances alone would show it never getting
    ## there, and once did.
    cases = [('gear2, ywr (default)', Gear2Integrator, 'ywr', 1e-3),
             ('gear2, ywr (default)', Gear2Integrator, 'ywr', 1e-4),
             ('gear2, classic', Gear2Integrator, 'classic', 1e-3),
             ('trapezoidal, ywr', TrapezoidalIntegrator, 'ywr', 1e-3),
             ('trapezoidal, ywr', TrapezoidalIntegrator, 'ywr', 1e-4)]

    print(".. list-table:: Stiff RLC from a charged initial state: accepted step"
          " ratios against the bound %.6f" % ZERO_STABILITY_RATIO)
    print("   :header-rows: 1")
    print("   :widths: 22 10 10 14 14 16")
    print("")
    print("   * - configuration")
    print("     - reltol")
    print("     - steps")
    print("     - force-accepts")
    print("     - worst ratio")
    print("     - ratios over bound")
    for label, cls, lte, reltol in cases:
        n, forced, worst, over = run(cls(lte), reltol)
        print("   * - %s" % label)
        print("     - %.0e" % reltol)
        print("     - %d" % n)
        print("     - %d" % forced)
        print("     - %.4f" % worst)
        print("     - %d" % over)

The last column is the gate: it must be zero on every row.  The rows that reach
the escape hatch at all are the ones to read -- before the repair, the default at
reltol 1e-3 produced one accepted step ratio of exactly 10.0, and the trapezoidal
row at 1e-4 produced eight above :math:`1+\sqrt 2`, the largest 6.4669.

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
