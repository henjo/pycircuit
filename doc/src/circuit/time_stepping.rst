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

Performance notes: where the time actually goes
================================================

This section exists so the next person does not re-derive it. The measurements are
from the 139-unknown leapfrog benchmark
(``benchmarks/transient_stage2.py``); reproduce with ``--reps N`` rather than
trusting the numbers, and see the caution about timing on a loaded machine below.

**Assembly dominates, not the linear solve.** Before the stage-2 work, building
:math:`i`, :math:`q`, :math:`u`, :math:`G` and :math:`C` was about **96%** of a
transient's runtime and the linear solve about **2%**. This is the opposite of the
intuition carried over from large sparse simulators, and it has a consequence worth
stating plainly: **optimising the matrix solver first would have been wasted work.**
At :math:`n \sim 10^2` the solve is already negligible.

Three changes, each behaviour-preserving:

**Hoisting the per-element probes.** ``hasattr(self.toolkit, 'add_at')`` was
evaluated once per element per stamp. ``NumericToolkit`` has no ``add_at``, so each
of those went through ``Toolkit.__getattr__``, which *formats an error message and
raises* — hundreds of thousands of exceptions per run to answer a question whose
answer cannot change during a solve. Attribute-lookup machinery was about **34%** of
total runtime; it is now under **5%**.

**Scattering once instead of per element.** ``np.add.at`` is numpy's unbuffered
scatter and is slow by design (it exists to handle duplicate indices correctly).
Calling it once per element pays that cost per element. Collecting the indices and
values and accumulating once with ``np.bincount`` is the *same sum in the same
order* — both walk the index array in sequence and add duplicates as they are met —
so the result is bit-identical, not merely equal to rounding. ``bincount`` is
real-only, so complex and object dtypes keep the ``np.add.at`` path.

**Not recomputing the charge vector.** The step controller and the history roll each
called ``cir.q(x)`` at a state the assembly had *just* evaluated. Measured 5.08
:math:`q` assemblies per accepted step against 3.06 for every other stamp; the
difference was exactly those two redundant calls. They are now served from a
value cached against the state that produced it.

**BLAS threads.** Circuit matrices are small enough that a threaded LAPACK spends
more time on thread overhead than on arithmetic — and the cost is not confined to
the solve, it is spread across the many small array operations in assembly. Limiting
BLAS to one thread measured **1.72x** on the whole transient. pycircuit will do this
itself if :mod:`threadpoolctl` is importable (see
:func:`~pycircuit.circuit.transient.blas_single_thread_available`); it is an optional
dependency, and without it the same win is available from the environment::

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python your_script.py

.. warning::

   **Timing on this class of machine is unreliable at the single-sample level.**
   Individual tests have been observed to vary by up to :math:`\pm 57\%` between
   runs, and a suite comparison once read as a 35% regression that controlled
   measurement showed to be 3%. Take the **minimum of several repetitions**, and
   check any claimed speedup against something that does not depend on the clock —
   step count, assembly count, or profile call count — before believing it.

.. note::

   One consequence of the BLAS thread limit: a threaded LAPACK reduces in a
   different order, so results differ in the last bits between one and many threads
   (measured: identical step count, :math:`\max|\Delta v| = 2.16 \times 10^{-15}` V).
   That is floating-point non-associativity in the BLAS, not a change of method, but
   it does mean a run is only bit-reproducible against another run with the same
   thread count.
