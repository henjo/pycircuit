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

Which devices take part
=======================

``Diode`` and ``BJT`` (stage 14). A device opts in by declaring
``pcnr_junctions``, a list of ``(anode, cathode)`` terminal pairs. Each pair
becomes one of the device's own unknowns, and the device is evaluated entirely at
them: its current, its Jacobian and its charge.

A BJT has two such unknowns, ``v_be`` and ``v_bc``. Its currents depend on both, so
PCNR evaluates a device once per iteration from the vector of its unknowns, and
each one's column in the enlarged Jacobian has an entry on every terminal (three,
not two). Any ``Semiconductor`` gets this from its own ``eval_i_pure`` and
``eval_q_pure``, provided its junction pairs form a spanning tree over its
terminals. The terminal voltages are then rebuilt from the unknowns, and the model
is still written once. A declaration that is not a spanning tree is refused when
the junctions are discovered.

**Charge.** In transient, a participating device's charge is evaluated where its
current is, at its own unknowns. Otherwise the charge would sit at the node
voltage while the current sits at ``v_lim``, which is the inconsistency PCNR
exists to remove. The rest of the circuit's charge is assembled without the
device, the device's ``q(v_lim)`` is added back before the integrator sees it,
and its ``alpha0 * dq/dv`` joins its Jacobian column (``geq = alpha0 * C`` for
every integrator). On a BJT switching stage, PCNR and classic limiting take the
same 1136 accepted and 36 rejected steps. The Jacobian given to the step
controller matches a finite difference of the full transient residual to 7.5e-10.

**What does not take part, and what that means.** ``JFET`` and ``Varactor`` have
classic limiters but no PCNR junctions. Under ``pcnr=True`` they are not limited
at all, and the same holds for any device without ``pcnr_junctions``. The
transient's own DC operating point still uses classic limiting.

What PCNR removes, stated correctly
===================================

**This page previously claimed that ``Diode.G`` linearises around the limited
voltage while ``Diode.i`` evaluates at the node voltage, so that "the Jacobian and
the residual are taken at different points". That is false**, and it mattered:
it made PCNR look like a correctness fix when it is not.

``Diode.i`` ends with ``I_eff = I - g * (VD - v_nodes)`` — the companion form
``i(x) = I(VD) + g(VD)(v(x) − VD)``, which is affine in ``x`` with slope exactly
``g(VD)``; and ``G`` returns ``g(VD)``. So ``G`` *is* ``di/dx``, measured by
central difference to a worst relative disagreement of 5.045e-10 across node
voltages from −2 V to +0.7 V. Re-measure with
``benchmarks/transient_review/companion_consistency.py``.

The residual and the Jacobian are consistent **within one linear solve**. What
PCNR actually removes is the *hidden* status of the linearisation point ``VD``:
it depends on the previous iterate and on the order devices are visited, and two
devices sharing a branch can undo each other's adjustment. Across iterations
``VD`` moves, so limiting is a modified Newton on a sequence of tangent points
rather than Newton's method on the original system.

**PCNR does not change the equations; it changes what is allowed to be
implicit** — the tangent point is promoted from a per-device side effect to an
unknown the solver carries and converges. Its benefit is convergence robustness,
not a different answer. ``Semiconductor`` addresses the same hidden-state problem
differently, by moving the evaluation point so everything is evaluated
consistently.

A companion far from the operating point is *valid*, not broken: at −2 V the
tangent above presents 6.77e-04 S where the device has 9.49e-47 S. That is why a
**stale** ``VD`` is silently a plausible matrix, and how the defect described
below survived three gates.

What it costs, and what it does not buy
=======================================

PCNR converges everywhere classic limiting does, including six circuits chosen to
make limiting matter, and agrees with it to 1e-6 relative. Iteration counts are
within one throughout, so the extra consistency buys no extra iterations.

The PCNR devices are left out of the ordinary assembly (``SubCircuit.G`` and
``SubCircuit.i`` take a ``skip`` set) and stamped once, at their own unknowns. On
60 diodes a PCNR Newton iteration costs 772 µs against 1081 µs for classic
limiting. Until gate 13-8 it cost 1774 µs: every PCNR device was assembled at the
node voltage and then subtracted out again. That only stayed finite because of the
trap described below. The JAX toolkit, which evaluates each device class in one
batch, honours ``skip`` by slicing the batch (gate 13-11). PCNR therefore runs on
a JAX-toolkit circuit with the numeric toolkit's iteration counts and answers.
Before 13-8 it went nan there on any large overshoot, because autodiff makes the
diode honest and the old subtraction then cancelled an overflowed current.

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

The trap that was here, and how it was closed
=============================================

``Diode`` keeps the limiter's output in ``_vlim`` and linearises around it. Only
``Diode.limit`` moves ``_vlim``, and only the three Newton loops that call
``cir.limit`` call that. Until gate 13-9, every *other* caller of ``cir.G`` or
``cir.i`` got the companion around whatever ``_vlim`` was last written. That is a
valid tangent of a frozen linearisation, with the type and magnitude of a
Jacobian. It surfaced twice:

- **PCNR's step controller**: ``cir.G(x)`` carried no diode conductance, which
  changed the transient step count by 6.6x (gate 13-6). Inside
  ``augmented_system`` the same wrong value was added and subtracted, so every
  DC test passed.
- **PSS**: its Newton has no limiter, so the diode was frozen at 0 V for the
  whole solve. A half-wave rectifier whose steady state peaks at 4.25 V came out
  as **0.000 V everywhere**.

Now the three limiting loops (``DC``, the standard transient and the classic
branch of the coupled path) run inside a classic-limiting scope. Inside it
``Diode`` behaves exactly as before. Outside it ``Diode`` is a plain function of
the node voltage and writes no state, so ``cir.G(x)`` is an honest Jacobian for
any caller. Classic DC and transient results are unchanged, bit for bit on all
but one test circuit. On a mains rectifier the waveform moves by 3.2e-14 with
the same steps: there, the step controller's Jacobian used to be taken up to
1.3e-12 V away from the converged point, and is now taken at it.

The order mattered. Making ``Diode`` honest **before** PCNR stopped
assembling-and-subtracting breaks PCNR outright. On a 1 A source the first
predict puts the node at 2.6e11 V, an honest diode overflows ``exp``, and the
subtraction is nan.

The coupled (Fang) path honours ``pcnr=True``. It solves for the step size
*inside* the Newton loop, so unlike the standard path it does **not** reproduce
classic limiting step for step. The LTE is evaluated at a mid-iteration iterate,
which is the limited point under limiting and the full update under PCNR, so the
two grids diverge slightly and legitimately. What must hold, and is measured, is
that neither is less accurate against an independent reference.
