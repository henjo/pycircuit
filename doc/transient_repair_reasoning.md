# Repairing the transient analysis — the reasoning

**Status: written 2026-07-30, before any stage has run. No outcomes recorded yet.**

This document says *why*. `transient_repair_plan.md` says *what*, in what order, and
what result would make each stage a success. Read this one first if you want to argue
with the decisions; read that one if you want to execute them.

## How we got here

None of this was the objective. The objective was to revive
`benchmarks/nonlinear_leapfrog_sweep.py`, a transient-vs-perturbation cross-check on
the 127-unknown leapfrog that had been committed but had never completed a run. The
recorded reason it never ran was a **cost wall** — "2 ms of simulated time did not
complete in 6 minutes" — attributed to stiffness.

Raising the drive from 1 kHz to 100 kHz should have fixed that by shortening `tend`
100x. It did not, and measuring why turned up four separate defects in the transient
engine, two of which are already fixed. The cost wall was real but it was not
stiffness; and behind it sat a problem that is not about cost at all.

**The lesson worth keeping, because it shaped everything below:** the recorded
diagnosis ("stiff, therefore slow") was plausible, was never measured, and was wrong.
The mean step was 5.03 ns against a 39 ns cap — a ratio that says *something is
rejecting steps*, not *something is stiff*. One measurement of the step-size
distribution would have found this at the time. Stages below therefore each declare
what they will measure, not only what they will change.

## What is broken, and how we know

**(A) `Gear2Integrator.compute_lte`, the `'classic'` branch, estimates the wrong
derivative.** `integrator.py:192-206`. BDF-2's truncation error needs `q'''` scaled by
`h**2`; the code takes a second divided difference of charge — which is `q''` — scaled
by `h**3`. That is not a rescaling slip: dimensionally it is not a current.

Evidence, at formula level with full healthy history so no first-step effect can
contribute: on `q(t)=sin(2pi*1e6*t)` the estimate is ~1e-15 of the true error and
scales as `h**2.99` instead of `h**2.00`. The consequence is not inaccuracy, it is a
**dead controller**: with `err ~ 3e-11` against an accept threshold of 1,
`min(2.0, 0.9*(1/err)**(1/3))` in `stepcontroller.py:79` saturates at 2.0 every step,
so `h` doubles to `max_step` and pins there with zero rejections, permanently. Results
are *bit-identical* across a 1e3 change in `abstol` and `reltol`. Reacting at all
would need tolerances ~3e10x tighter.

Why the code does this at all is understandable: Gear-2 keeps only two past charges
(`get_required_history() == 2`), so a third divided difference of `q` is unavailable.
The `'ywr'` branch sidesteps it using the companion-current history. `'classic'`
silently dropped a derivative instead.

**This is already documented.** `doc/src/circuit/lte_dae.rst:154-165` says the classic
estimate "under-estimates the true truncation error" and that "the controller is then
fooled into taking too few steps". Two things are nonetheless wrong: `Gear2Integrator()
still defaults to it`, which is why it was hit at all; and the docs frame the impact as
"roughly an order of magnitude" of accuracy when the controller is in fact entirely
disabled. A known defect left as the default is a worse failure mode than an unknown
one, because the documentation reads as though someone decided it was acceptable.

> **Fixed 2026-07-30.** Stage 1 replaced the estimate; stage 2b moved
> `Gear2Integrator`'s default to `'ywr'`, so the sentence above about the default
> no longer describes the tree. Stage 5 corrected the `lte_dae.rst` claim.
> Measurements: `transient_repair_plan.md`.

**(B) A test asserts the bug.** `test_analysis_transient.py::test_lte_formula_ywr`
asserts `n_gy > n_gc` — that `'ywr'` takes more steps than `'classic'`. That is true
only because `'classic'` takes the fewest steps physically possible. With a correct
estimate the assertion fails. So (A) cannot be fixed without a decision about (B),
which is why (A) has not been applied.

**(C) The same expression is copied into the JAX path.** `jaxtransient.py:247-254`
(`estimate_lte`, `method='gear'`) is line-for-line the same defect, and its Trapezoidal
branch (257-260) has the same missing-derivative shape. Dormant, because that path
defaults to `'ywr'` — but dormant is not fixed, and the copy means a future reader
finds two mutually-confirming wrong implementations.

**(D) Post-breakpoint steps are accepted with no error check, and that is periodic.**
`IntegralController.evaluate_step` returns `(True, h_curr)` unconditionally when
`is_first_step` (`stepcontroller.py:26-28`; it also assigns an `err = 0.5` that is
never used). `transient.py:314` re-arms `_is_first_step` after **every breakpoint**,
and `Sin.next_event` (`func.py:30`) fires every quarter period. So a `VSin` drive
produces a periodic, drive-synchronous, full-`max_step` step that is never checked.

This affects **all** integrators, not just Gear2, and it converts `max_step` from a
cost knob into a correctness knob — which is not what a user reading `timestep=` would
expect. It is also the reason the original 2-step Gear2 observation was uninformative
on its own: step 1 was the unconditional accept and step 2 was the only checked step
in the run.

**(E) Nothing in the test suite could have caught any of this.** Twelve tests touch
Gear2. `test_vss_gear2.py` never calls `compute_lte`. `test_transient_methods_step_response`
uses `fixed_timestep=True`, so the LTE path is not entered.
`test_transient_adaptive_efficiency` asserts `adapt_steps < fixed_steps` — which a
*blind* controller satisfies most emphatically of all, by taking the fewest possible
steps. Nothing anywhere asserts that a step is ever rejected, that step count or
accuracy is monotone in `reltol`/`abstol`, or that `compute_lte` scales with the right
power of `h`. The third of those is a fast unit test needing no circuit at all.

**(F) Already fixed, for context.** `e37ddad`: the step controller was handed the
residual-flavoured tolerance vector (`iabstol` on node rows) but applies it to
`lte = J^-1 Eg`, which is in solution units — so every node row used 1 pA as a voltage
tolerance; and `vabstol` defaulted to 1e-12 V against Spectre's and SPICE VNTOL's 1 uV.
19x fewer steps, full suite at the 715-test baseline.

## Scope

### In

1. **(A) Fix `compute_lte`.** A verified replacement exists: second divided difference
   of `g` from the `iq` history, asymptotically exact (estimate/true ratio 1.0115 ->
   1.0007 as `h` halves, against `'ywr'`'s constant 0.75 bias). This is the defect that
   disables the controller, so it is the point of the exercise.
2. **(B) Rewrite `test_lte_formula_ywr`** to assert a property that is true of a
   *correct* implementation rather than of the current one.
3. **(D) Fix the unconditional first-step accept.** Narrowly: the genuine first step of
   a run has no history and cannot have its LTE estimated, so accepting it is correct;
   re-arming that exemption at every breakpoint is not. The fix is to distinguish "no
   history" from "history was just reset", and to bound the step in the latter case.
4. **(E) Add the three missing test classes**: `compute_lte` scales with the right
   power of `h` (unit, no circuit); steps are actually rejected on a stiff case;
   step count and accuracy are monotone in `reltol`.
5. **(C) Fix the JAX copy** of the same expression, or delete the dead branch.
6. **Correct `lte_dae.rst`** to say the controller is disabled rather than an order of
   magnitude less accurate.

### Out

- **The leapfrog instability.** `leapfrog_5th_order` appears to have two right-half-plane
  poles (s = +1.4491e+05, +5.6716e+04). That is a *fixture* defect, not a transient-engine
  defect, and it is under separate investigation. It is out of scope here because the
  engine must be trustworthy before it can be used to judge a circuit — fixing them
  together would leave neither verified.
  **Reconsider if** the instability turns out to be caused by the engine (e.g. an
  integrator sign error manifesting as a growing mode) rather than by the topology.
- **Changing the `reltol` default from 1e-4 to Spectre's 1e-3.** Defensible either way,
  and unlike `vabstol`'s 1e-12 it is not indefensible, so it is a preference not a bug.
  **Reconsider if** stage 4's monotonicity tests show 1e-4 is being routinely
  overridden, or if suite runtime becomes the binding constraint.
- **PSS / shooting (`shooting.py`) instead of long transients.** The right tool for a
  periodic steady state in principle. Rejected for now because a shooting Newton needs
  the monodromy matrix, and if that is finite-differenced over 127 unknowns it costs
  ~127 single-period integrations per iteration — worse than the ~20 periods a plain
  transient needs.
  **Reconsider if** `PSS` turns out to use a matrix-free or directly-computed
  sensitivity, in which case it is strictly better and the whole cost analysis changes.
- **The JAX transient as a speed lever.** `jaxtransient.py` has jit-compiled Newton
  loops and could plausibly cut the ~0.2 s/step Python MNA assembly by 10-100x, which
  would make the leapfrog comparison routine rather than overnight. Out of scope
  because it is a performance project, not a correctness one, and correctness first.
  **Reconsider if** the repaired engine is correct but the leapfrog comparison is still
  unaffordable — at which point this becomes the only remaining lever.
- **Reworking `TRTOL = 7.0` and `MAX_REJECT = 3`**, which are hard-coded rather than
  parameters (`transient.py:260`, `:269`). Noted, not touched: changing accept
  thresholds while also changing the error estimate would make both unattributable.
  **Reconsider if** stage 1 lands and the rejection rate is then pathological.

## Reviewed as an architect: what does this force to change outside itself?

The lens matters here because stage 1 is a one-function change with a very wide shadow.

**Every step count in the project changes.** A controller that currently never rejects
will start rejecting. Any test asserting a step count, a wall-clock time, or an
"adaptive beats fixed" inequality is at risk — `test_transient_adaptive_efficiency` is
explicitly built on such an inequality, and `test_lte_formula_ywr` is already known to
fail. So stage 1 must be run against the *whole* suite, not the transient tests, and
the plan cannot assume the blast radius is one file.

**The suite's runtime is already dominated by one transient test.**
`test_stress_stiff_rlc_pulse` is ~266 s of a ~492 s total (recorded as P15 in
`architecture.md`). If a correct LTE estimate makes Gear2 take more steps, this can get
worse, and a repair that doubles suite runtime will be resented into being reverted.
This is the main argument for measuring runtime as a declared gate rather than
noticing it afterwards.

**`Gear2Integrator`'s default is a public choice, not an implementation detail.**
Switching the default to `'ywr'` fixes the symptom in one line and is tempting.
It is the wrong primary move: it leaves a wrong `'classic'` in the tree for the next
reader, and `'ywr'` carries a constant 0.75 bias that a correct `'classic'` does not.
Fix the formula; then decide the default on its merits.

**Two implementations must not drift.** `jaxtransient.py` duplicates the LTE logic
rather than sharing it. Fixing one and not the other is how this defect survived in
two places to begin with, so (C) is in scope even though that path is dormant.

## What could make this whole thing wrong

- The verified fix for (A) was checked against an analytic `q(t)` and on small stiff
  circuits, not on the 127-unknown leapfrog. Agreement there is an assumption, and
  stage 1's gate must not be stated in terms of the leapfrog.
- The 40% waveform discrepancy originally attributed to Gear2 on the leapfrog is
  probably compounded by (D) rather than purely a Gear2 effect. If stages 1 and 3 are
  done together, the attribution will not be recoverable — hence they are separate
  stages with separate measurements.
- Fixing (D) changes behaviour for every source with breakpoints, including `VPulse`,
  where the breakpoint exists precisely because there *is* a discontinuity and a small
  step is genuinely wanted. The fix must not turn a legitimate discontinuity into a
  rejection loop.
