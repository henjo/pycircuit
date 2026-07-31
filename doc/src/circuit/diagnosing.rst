===========================
Diagnosing a failed run
===========================

A simulation can fail in three ways, and they want different responses.  This page
is about telling them apart from what the simulator says.

Before stage 6 it could not tell you: a floating node, a genuinely non-convergent
bias point and an ill-conditioned matrix all produced the same sentence --
``Source Stepping failed at lambda=0.0`` -- which names the last thing that was
tried rather than the first thing that went wrong.

The structure is not incidental
===============================

Continuation methods -- gmin-stepping and source-stepping -- exist to guide Newton
to a solution that *exists* but is hard to reach.  They cannot conjure an equation
that is missing.  So the simulator classifies a **structurally** singular system
before continuation is attempted, and reports it as
:class:`~pycircuit.circuit.analysis.SingularMatrix` rather than as a convergence
failure.  Those two exception types are the first thing to read.

.. list-table:: What the exception type already tells you
   :header-rows: 1
   :widths: 26 34 40

   * - Exception
     - Meaning
     - What to do
   * - ``SingularMatrix``
     - The system has no unique solution *as written*.
     - Fix the circuit. Continuation will not help and is not tried.
   * - ``NoConvergenceError``
     - A solution exists; Newton could not reach it.
     - Look at the named node, then at tolerances and the initial guess.

"No DC path to ground"
======================

The commonest structural failure.  Every node needs a conductive path to the
reference node at DC -- a node reachable only through capacitors is not
determined by the equations, however sensible the circuit looks.

.. code-block:: text

    SingularMatrix: singular Jacobian: 'floaty' appears in no equation, so
    nothing determines it -- for a node that means no DC path to ground
    (add a resistor, or use uic=True to skip the operating point)

The node is named.  The two fixes are a large resistor to ground, or ``uic=True``
to skip the operating point entirely and start the transient from zero -- which is
the right answer when the circuit *has* no DC solution by design.

A convergence failure names the worst unknown
=============================================

.. code-block:: text

    Source Stepping failed at lambda=1.0: Gmin Stepping failed at gmin=0.001:
    StandardNewton failed to converge after 2 iterations;
    residual worst at 'c': |f| = 0.01709 against a tolerance of 1.709e-06 (1e+04x over);
    update worst at 'c': |dx| = 1.907 against a tolerance of 5.933e-05 (3.21e+04x over)

Read it right to left.  The rungs on the left are the continuation ladder --
useful context, and nothing more.  The part after the last colon is the diagnosis:
**which** unknown failed, on **which** of the two convergence tests, and by how
much relative to its own tolerance.

Both tests must pass.  ``|f|`` is the KCL residual at that node -- current that
does not balance.  ``|dx|`` is how much the solution was still moving.  A large
``|dx|`` with a small ``|f|`` usually means a nearly flat region; the reverse
usually means the step is being limited.

The multiplier is the useful number.  ``1e+04x over`` says the run is nowhere
near; ``1.4x over`` says a slightly looser ``reltol`` or a few more ``maxiter``
would probably do it.

What the run actually did
=========================

A completed transient carries a statistics object, on both the analysis and the
result:

.. code-block:: text

    >>> res = Transient(cir).solve(tend=3e-6, timestep=1e-7)
    >>> print(res.statistics)
    accepted 153, rejected 18 (10.5% of attempts), Newton iterations 331 (2.2 per accepted step)
    force-accepts 0, order drops 20, breakpoints hit 26
    step 1e-10 .. 6.009e-08 s
    time 0.225 s total, 0.196 s in the Newton solve (87.2%)

Use ``res.statistics.as_dict()`` for the raw values.

**Read ``force-accepts`` first.**  It counts steps the solver accepted whose
truncation error it could not bound -- it had rejected the step three times, and
took it anyway so that time would advance.  Those steps also emit a
``RuntimeWarning`` naming the time.  On every circuit measured in this project the
count is now zero, so **a non-zero value is the run telling you that part of its
own result is not error-controlled.**

The rest read as you would expect, with two worth calling out:

* a **rejection rate** above roughly a third means the step controller is
  fighting the estimator rather than following it;
* **min step** orders of magnitude below max step is normal near a discontinuity
  and suspicious if the drive is smooth.

Order drops are not a problem in themselves: the integrator deliberately drops to
first order across a breakpoint, so a ``VPulse`` drive produces one per edge.
