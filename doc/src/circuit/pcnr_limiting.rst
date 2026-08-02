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

What PCNR would change, and why it is not this
==============================================

The paper — Aadithya, Keiter & Mei, *"Predictor/Corrector Newton-Raphson
(PCNR)"*, Sandia — makes three changes, none of which is present here:

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

The inconsistency PCNR removes is real and is present here. On a toolkit without
automatic differentiation, ``Diode.G`` linearises around the *limited* voltage
while ``Diode.i`` evaluates at the *node* voltage — the Jacobian and the residual
are taken at different points. ``Semiconductor`` avoids it differently, by moving
the evaluation point itself so that everything is evaluated consistently.

Whether PCNR is adopted depends on two measurements that have not been made:
whether the Schur elimination is genuinely free, and whether the method converges
at least as often as limiting does on the nonlinear stress circuits. **A method
that is more consistent and converges less often is not an improvement.**
