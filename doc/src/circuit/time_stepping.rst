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

Option A: ``coupled_lte=True`` — what it actually is
----------------------------------------------------

.. warning::

   **This section previously described a method that is not implemented.** It
   claimed an augmented :math:`(N+1)` system, an analytical gradient
   :math:`E_h = p(E + TRTOL)/h`, and a "golden window"
   :math:`0.7\tau \le \epsilon \le 3.0\tau` following Fang §3.3. A review against
   the source paper in 2026-07 found **none of those three in the code**. The
   description below is what ``_solve_coupled`` does. The discrepancy is recorded
   rather than quietly deleted, because a documented-but-absent feature is the
   same class of defect as a silent wrong answer: both make a claim the software
   does not honour.

What Fang's paper proposes (DAC 2013, §3.1–3.4) is genuinely a co-determination
method: the step size :math:`h_m` is an **unknown**, solved together with the
circuit equations as :math:`N+1` nonlinear equations via a bordered system

.. math::
    \begin{pmatrix} J & p \\ q^T & d \end{pmatrix}
    \begin{pmatrix} \Delta v \\ \Delta h \end{pmatrix} =
    \begin{pmatrix} -f_{ckt} \\ -f_{lte} \end{pmatrix}

with :math:`p = \partial f_{ckt}/\partial h_m`, :math:`q^T = \partial f_{lte}/\partial v_m`
and :math:`d = \partial f_{lte}/\partial h_m`. Its §3.4 "approximate Newton"
variant avoids re-solving by *correcting* the solution already computed,
:math:`\Delta v^{k+1} = \Delta v^{k+1/2} - J^{-1} p (h^{k+1} - h^{k})`.

``_solve_coupled`` implements neither. It converges the circuit at :math:`h`,
evaluates the LTE, and if the LTE is over tolerance it shrinks :math:`h` and
**re-solves from scratch**, up to ``MAX_LTE_ITERS = 10`` times. That is a
rejection loop — structurally the same scheme as Option B, with a different retry
limit. There is no :math:`p`, no bordered system, and no error window; the
``analytical_eh`` argument that survives in the signature is a vestige of the
:math:`E_h` gradient described above, and is never read.

Consequently ``coupled_lte=True`` does **not** eliminate rejected steps, and the
two options are not the algorithmic contrast this page once claimed. Prefer the
default (``coupled_lte=False``).


Starting point: what happens when the operating point fails
============================================================

A transient needs an initial state. Unless you supply one, it is the DC operating
point, solved by an inner :class:`~pycircuit.circuit.dcanalysis.DC` constructed
from the transient's own toolkit, environment parameters, tolerances, solver and
scaler — so the bias point is found under the same conditions as every step that
follows it.

**If that solve fails, the transient raises.** It does not substitute zeros.

This matters more than it sounds. Until 2026-07 a failed operating point was
silently replaced by a vector of zeros, and the run continued to completion:
a circuit that had never been biased produced a full, plausible-looking waveform,
and nothing in the result distinguished it from a correct one. That is the most
expensive kind of defect a simulator can have — it does not fail, it lies.

Two circuits legitimately have no operating point, and both are common:

* an **ideal integrator** (``Idt``, ``Idtmod``) — its output is the unbounded
  integral of a constant input, so no steady state exists;
* a **charge pump** or any topology whose nodes reach ground only through
  reverse-biased junctions or capacitors — structurally singular at DC.

For these, say so:

.. code-block:: python

    # start from zeros, deliberately (SPICE's "use initial conditions")
    tran = Transient(cir, uic=True)

    # or start from a state you computed yourself
    res = tran.solve(x0=my_operating_point, tend=..., timestep=...)

Both are explicit choices. The error raised when the DC fails names them.

Tolerances: which knob does what
=================================

Two distinct roles, and until 2026-07 they shared one parameter:

``vabstol``, ``iabstol``, ``reltol``
    **Newton's** convergence criteria, shared with :class:`DC` so the operating
    point and the steps after it are solved to the same accuracy.

``lte_vabstol``, ``lte_iabstol``
    The **step controller's** tolerances, applied to the local truncation error
    :math:`J^{-1}E_g`, which carries the units of the solution vector.

Changing ``lte_vabstol`` moves the step count; changing ``vabstol`` does not.
They were one parameter, which meant relaxing the step controller silently
loosened Newton's node convergence by a factor of :math:`10^6` — while ``DC``
kept the tighter value, so the operating point was solved a million times more
precisely than any step that followed it.

.. note::

   ``lte_vabstol`` is an interim measure. Spectre carries a single tolerance set
   and derives the LTE bound by multiplying it by ``lteratio``; what pycircuit is
   missing is ``relref``, the choice of *reference* for the relative term. Under
   pycircuit's fixed per-node reference, a node carrying no signal has its
   tolerance collapse to the absolute floor, which is what forced the absolute
   floor upward in the first place. See ``doc/transient_work_plan.md``.
