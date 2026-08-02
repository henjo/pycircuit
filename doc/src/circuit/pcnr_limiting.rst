===============================================================
Newton limiting for exponential devices
===============================================================

.. note::

   **This page previously claimed that pycircuit implements PCNR. It does not.**
   The claim was reviewed against the paper and the source on 2026-08-02 and
   found false in all three of its specifics: there are no extra unknowns for
   limited quantities, there is no prediction/correction split, and the only
   Schur-complement code in the tree is ``SchurCoupledNewton``, which solves for
   the *time step* (Fang, DAC 2013) and has nothing to do with limiting.

   The page now describes what is implemented. PCNR remains a candidate and is
   scheduled as STAGE 13 in ``doc/transient_work_plan.md``; the section below
   says what adopting it would change. A documented-but-absent feature is the
   same class of defect as a silent wrong answer, and this is the second one
   found in this project — the first was the Fang citation in ``_solve_coupled``.

The problem
===========

Newton-Raphson linearises, and an exponential device punishes that. A step taken
from a lightly-biased point overshoots enormously, and the model is then
evaluated somewhere it has no business being. Clamping the exponent stops the
result becoming ``nan``; it does not stop the iteration wandering. Measured: with
the clamp alone and no limiting, a common-emitter stage converges to a **genuine
but spurious** operating point 200 V below the rails, with the base-collector
junction forward biased at 0.91 V carrying the 20 A the base resistor delivers.
That is a solution of the circuit equations; it is simply not the one anyone
wants.

What pycircuit implements
=========================

Classic SPICE ``pnjlim``, bounding the per-iteration excursion of a junction
voltage. ``vc`` is the voltage at which the exponential's curvature starts to
dominate; above it the update is compressed logarithmically. **The limiter only
moves the point the next Jacobian is taken at — it never changes the equations**,
so the converged answer is a solution of the unmodified system.

Devices declare their junctions rather than implementing a limiter::

    junctions = ((anode, cathode, move), ...)   # indices into `terminals`

``move`` is not redundant. A BJT's two junctions share the base as anode, so
limiting both by adjusting the anode would have the second undo the first;
moving the base for B-E and the collector for B-C keeps them independent.

An empty ``junctions`` means no limiting, and for ``ZenerDiode`` that is
deliberate rather than an omission: a forward-junction limiter would step
straight through the breakdown knee, which is a second exponential running the
other way. It gets the direction-agnostic exponent clamp instead.

What PCNR changes, and what it is worth
=======================================

The paper — Aadithya, Keiter & Mei, *"Predictor/Corrector Newton-Raphson
(PCNR)"*, Sandia — makes three changes. All three are now implemented in
``pycircuit/circuit/pcnr.py``, selectable with ``pcnr=True`` and **off by
default**:

1. **Every limited quantity becomes an MNA unknown.** The vector grows to
   ``x = [x_MNA; x_lim]``, with a residual row ``v_Dk − (e_a − e_b) = 0`` tying
   each new unknown to the branch voltage it stands for. Devices then evaluate
   at *their own* ``v_Dk``.
2. **Each NR iteration splits into predict and correct.** Predict is a full
   Newton step on the enlarged system; correct is each device limiting **only the
   variables it owns**, which is why devices cannot interfere with one another.
3. **The extra unknowns are eliminated by a Schur complement**, which is cheap
   because the ``lim/lim`` block is the identity, so the linear system actually
   solved stays the size of the original MNA system.

The inconsistency PCNR removes is real. On a toolkit without automatic
differentiation, ``Diode.G`` linearises around the *limited* voltage while
``Diode.i`` evaluates at the *node* voltage — the Jacobian and the residual are
taken at different points. ``Semiconductor`` avoids it differently, by moving the
evaluation point itself so that everything is evaluated consistently.

What it costs, and what it does not buy
=======================================

PCNR converges everywhere classic limiting does, including six circuits chosen to
make limiting matter, and agrees with it to 1e-6 relative. Iteration counts are
within one throughout, so the extra consistency buys no extra iterations. The
cost is per-iteration assembly work, and it is not free: devices still stamp
their nonlinear part into ``G`` and it has to be subtracted again. Reaching the
paper's cost would mean devices not stamping that part at all — an assembly
change, and the one this implementation deliberately avoided in order to remain a
layer.

**In transient, PCNR and classic limiting produce identical results** — the same
step count, the same step times, and the same waveform to nine figures on a
half-wave rectifier, a full-wave rectifier, a charge pump and a clipper.

That identity is not a disappointment; it is the correct answer, and worth stating
plainly because an earlier version of this page claimed otherwise. The limiter and
PCNR's ``refine`` change only the *iteration path*. A converged Newton solution is
whatever the residual says it is, independent of the route taken to it, so equal
solutions give equal truncation error and hence equal steps. **A measured
difference in transient step count between the two paths is a defect signature,
not a feature.** It is asserted as such by
``test_gate_13_6_pcnr_and_limiting_take_the_same_steps``.

A trap for anyone extending this
================================

``Diode.G`` linearises around ``self._vlim``, which only ``Diode.limit`` writes —
and PCNR never calls ``limit``, since limiting is what it replaces. During a PCNR
run ``_vlim`` therefore holds a stale value and **``cir.G(x)`` returns a Jacobian
with no diode conductance in it**.

Inside ``augmented_system`` this cancels exactly: the same wrong value is added by
``cir.G(x)`` and subtracted again by the per-device loop, so the system that gets
solved is correct and every converged voltage is right. That is why the defect
survived every DC test. It escaped through one door only — the matrix handed to
the transient step controller, which forms ``lte = J⁻¹ Eg`` — and there it changed
the step count by 6.6x.

So: **do not call ``cir.G(x)`` on the PCNR path and expect a usable Jacobian.**
The matrix to use is the Schur matrix that ``predict`` already factorises; at
convergence ``v_lim == e_a − e_b``, so it *is* the Jacobian of the residual with
respect to ``x``. The underlying hazard is only contained, not removed — the clean
fix is the assembly change described above.
