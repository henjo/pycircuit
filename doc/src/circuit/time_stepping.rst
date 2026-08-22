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

Option A: ``coupled_lte=True`` — Fang's coupled time-stepping
-------------------------------------------------------------

.. note::

   **This section described an absent feature for a long time, then described
   its absence.** It first claimed an augmented :math:`(N+1)` system, an
   analytical gradient and a "golden window" following Fang §3.3; a review in
   2026-07 found none of the three in the code and the page was rewritten to say
   so. As of stage 12B (2026-08-01) the method is implemented, and this section
   describes what runs. The history is kept because a documented-but-absent
   feature is the same class of defect as a silent wrong answer.

Fang's method (DAC 2013) makes the step size :math:`h_m` an **unknown**, solved
together with the circuit equations as :math:`N+1` nonlinear equations. Its
Figure 3 has no rejection branch at all: the predicted step is only an initial
guess, and *"there is no backup due to LTE"*.

The local truncation error is eq (6),

.. math::
    \epsilon_m = \left| v_i(t_m) - v_{i,\mathrm{extrapolated}} \right|

— a **solution-space** quantity: the computed solution minus a polynomial
extrapolation from previous time points, maximised over the unknowns, with
:math:`i` the *controlling LTE node*. This is **not** the charge divided
difference the standard path uses. The distinction is the whole reason the
method works: a charge divided difference carries repeated :math:`1/h` factors,
so it *grows* as the step shrinks and cannot be solved for :math:`h` — Newton
runs away from the root. See ``pycircuit/circuit/stepcontroller.py``'s
``SolutionLTEController`` and ``doc/fang_stage12_conclusions.md`` §3.

``Transient.fang_timestep`` implements Figure 4's two-stage Newton: solve the
ordinary :math:`N` system, estimate the LTE, and adjust the step **only if the
LTE condition fails**. The step correction is §3.4's approximate Newton, eqs (17)
and (18), which takes the new step from the error ratio and corrects the solution
already computed using the factors from the first stage.

**The bordered system of eq (12) is deliberately not used.** Its eq (14)
denominator :math:`q^T \Delta v_h + d` is the solution's sensitivity to the step
size minus the extrapolation's slope; both are approximately :math:`dv/dt`, so
the difference is the truncation error's derivative and is tiny by construction.
Measured at :math:`h = 1.6\times10^{-7}`: :math:`+1.818\times10^{9}` against
:math:`-1.820\times10^{9}`, a denominator of :math:`-2\times10^{6}` — three
digits lost, with the *sign* decided by the cancellation. This is very likely
what §3.4 means by the coupled system being "very sensitive to the change of step
size".

What it delivers, measured against closed-form solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Zero LTE rejections**, at every tolerance tried, on every circuit tried.
  That is Figure 3's central claim and it transfers intact.
* **Mean and median waveform error 10–20% lower** than the standard controller at
  matched ``reltol``. The *maximum* is a wash, so a comparison reporting only a
  maximum would miss this.
* **25–28% more Newton solves.** More accurate per unit tolerance, less accurate
  per unit work.

The paper's headline — 39% fewer time points — does **not** appear here, and the
reason is understood rather than mysterious: it comes from moving steps *up* onto
the tolerance, and ``IntegralController``'s prediction law is already deadbeat,
placing 91–96% of accepted steps within 5% of its target. The comparison method
in §4.1 redid a step only above a normalised LTE of 4.63 without resizing toward
it; pycircuit's baseline resizes.

Limitations
~~~~~~~~~~~

* ``TLine`` cannot be used with ``coupled_lte=True``. Its source vector comes
  from history interpolated at :math:`t - T_D`, so :math:`\partial u/\partial t`
  is the derivative of that interpolation and is unwritten; it raises rather than
  contributing a silent zero to :math:`p`.
* ``fixed_timestep`` and a caller-injected step controller are not honoured on
  this path. Breakpoints and ``uic`` are.
* No wall-clock comparison has been made; the figures above count Newton solves.

``coupled_lte=False`` remains the default. On this evidence it should be: the
coupled path buys zero rejections and slightly better average accuracy for
roughly a quarter more work.


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

* an **ideal integrator** (``Idt``, ``Idtmod``) *without an* ``ic`` — its
  output is the unbounded integral of a constant input, so no steady state
  exists.  Giving the element an ``ic`` pins its DC solution instead (the
  Verilog-A semantics; see :ref:`idtmod`) and removes the singularity;
* a **charge pump** or any topology whose nodes reach ground only through
  reverse-biased junctions or capacitors — structurally singular at DC.

For these, say so:

.. code-block:: python

    # start from zeros, deliberately (SPICE's "use initial conditions")
    tran = Transient(cir, uic=True)

    # or start from a state you computed yourself
    res = tran.solve(x0=my_operating_point, tend=..., timestep=...)

Both are explicit choices. The error raised when the DC fails names them.

The first step, and why ``reltol`` used to do nothing
======================================================

A run opens at ``firststep``. That step is the only one in the whole run whose error
**cannot be checked**: the truncation-error estimate is a difference against previous
points, and at the start there are none, so the controller has to accept it
unevaluated.

Until 2026-07 the opening step was ``timestep`` — which is also ``max_step``, the
largest step the controller is ever allowed to take. So the one step nobody could
check was also the biggest one, and its error dominated everything after it.

The symptom was that ``reltol`` did not control accuracy at all. On an RC step
response measured against its analytic solution, the global error was
:math:`1.3212 \times 10^{-1}` at ``reltol`` of 1e-3, 1e-4, 1e-5 **and** 1e-6 —
identical to five digits — while the step count went from 24 to 195. Eight times the
work for the same answer. Backward Euler and the two second-order methods likewise
agreed to five digits, which is why the choice of integrator appeared not to matter.

Opening at ``timestep * 1e-3`` instead costs one cheap step and lets the controller
grow from there. The same sweep now gives 3.13e-3 down to 3.44e-5, a **90.8x**
reduction, falling monotonically; and trapezoidal becomes 9.5x more accurate than
backward Euler while using 7.5x fewer steps.

The default is the knee of a measured curve, not a guess:

.. list-table::
   :header-rows: 1

   * - ``firststep``
     - steps
     - how much ``reltol`` can move the error
   * - ``timestep``
     - 50
     - **1.0x** — nothing
   * - ``timestep*1e-1``
     - 58
     - **1.0x** — still nothing
   * - ``timestep*1e-2``
     - 62
     - 63.9x
   * - ``timestep*1e-3`` (default)
     - 66
     - **90.8x**
   * - ``timestep*1e-6``
     - 75
     - 91.7x

Set it yourself if you need to::

    tran = Transient(cir, firststep=1e-12)   # explicit
    tran = Transient(cir)                    # None -> timestep*1e-3

.. note::

   The ramp is a **fixed** cost of roughly 15 steps — the geometric climb back to
   ``max_step`` — not a proportional one. On a long run that is negligible (+3.1% on
   the leapfrog benchmark, +4.2% over 300 time constants); on a run of only ten time
   constants it is +32%. Short runs pay proportionally more, and they are also the
   runs where the unchecked opening step did the most damage.

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
