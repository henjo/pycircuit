# Stage 12 — Fang's coupled time-stepping: what was found

**Status 2026-08-01.** The method is implemented and runs: `coupled_lte=True` solves for the
solution and the step size together and takes **zero LTE rejections**, which is the central
claim of the paper. It costs about 25% more Newton solves than the standard adaptive
controller on the circuits measured here. It is **not** the default.

Companion documents: `doc/fang_dac2013_math.md` (the paper's math, extracted from rendered
pages) and `doc/transient_work_plan.md` STAGE 12 (the plan and its gates).

Every number below is reproducible from `benchmarks/transient_review/stage12*.py`.

---

## 1. The headline: what the method delivers here, and what it does not

`benchmarks/transient_review/stage12b_coupled.py`, against closed forms (never a fine-mesh
reference — comparing the integrator with itself cannot detect a step count bought by
integrating less accurately). Opening step excluded; it is accepted unevaluated and no
tolerance governs it.

**rc-vsin** (RC driven by a sine):

| reltol | path | steps | reject | wall (s) | µs/point | mean err | median err | max err |
|---|---|---|---|---|---|---|---|---|
| 1e-5 | standard | 1288 | 19 | 2.750 | 2104 | 1.219e-03 | 1.423e-03 | 2.014e-03 |
| 1e-5 | **coupled** | 1591 | **0** | 2.828 | **1778** | **1.062e-03** | **1.034e-03** | 2.168e-03 |
| 1e-6 | standard | 4067 | 22 | 10.113 | 2473 | 3.830e-04 | 4.464e-04 | 6.350e-04 |
| 1e-6 | **coupled** | 5136 | **0** | 10.966 | **2135** | **3.298e-04** | **3.289e-04** | 6.626e-04 |

**stiff-rlc** (ringdown, `LC v'' + RC v' + v = 0`):

| reltol | path | steps | reject | wall (s) | µs/point | mean err | median err | max err |
|---|---|---|---|---|---|---|---|---|
| 1e-5 | standard | 490 | 4 | 1.168 | 2364 | 4.750e-03 | 3.869e-03 | 1.335e-02 |
| 1e-5 | **coupled** | 619 | **0** | 1.477 | 2386 | **4.134e-03** | **3.102e-03** | 1.335e-02 |
| 1e-6 | standard | 1499 | 5 | 3.237 | 2152 | 3.844e-03 | 2.948e-03 | 1.345e-02 |
| 1e-6 | **coupled** | 1954 | **0** | 4.238 | 2169 | **3.384e-03** | **2.419e-03** | 1.351e-02 |

**rc-pulse** (RC driven by a trapezoidal pulse — the only circuit here with real
breakpoints, and the one that found the defect in §5 item 9):

| reltol | path | steps | reject | wall (s) | µs/point | mean err | median err | max err |
|---|---|---|---|---|---|---|---|---|
| 1e-5 | standard | 1498 | 27 | 3.922 | 2572 | 9.577e-04 | 1.015e-03 | 2.063e-03 |
| 1e-5 | **coupled** | 2222 | **0** | 7.524 | 3386 | **6.607e-04** | **7.034e-04** | **1.414e-03** |
| 1e-6 | standard | 4537 | 29 | 9.311 | 2039 | 3.138e-04 | 3.303e-04 | 6.508e-04 |
| 1e-6 | **coupled** | 6584 | **0** | 15.804 | 2400 | **2.167e-04** | **2.297e-04** | **4.517e-04** |

**Conclusions.**

1. **Zero LTE rejections everywhere.** Figure 3 has no rejection branch and the
   implementation has none either. This is the claim of the paper that transfers intact.
2. **Better mean and median error at matched tolerance on all three circuits**, by 10–30%.
   On `rc-pulse` the *maximum* is better too, by 31%; on the two smooth circuits the maximum
   is a wash. Reporting only a maximum would have hidden this.
3. **It costs 23–45% more time points**, and **8–70% more wall clock**: +8% on `rc-vsin`,
   +31% on `stiff-rlc`, +70% on `rc-pulse`.
4. Net: **more accurate per unit tolerance, less accurate per unit work** — but the margin
   is circuit-dependent, and on the smooth sine drive the coupled path is within 8% of the
   standard path's runtime while eliminating its rejections.

> **Cost is measured per NEWTON ITERATION, not per time point.** An earlier version of this
> document divided wall clock by time points, because `newton_iterations` was flat zero on
> the coupled path — `fang_timestep` runs its own loop and never calls `_newton`, where the
> counter lives. That unit mixes two different things and cannot separate them.
>
> Measured properly: **the coupled path is 9–24% cheaper per Newton iteration on every
> circuit** (709 vs 781 µs on rc-vsin, 595 vs 786 on stiff-rlc, 694 vs 852 on rc-pulse). Its
> entire cost is doing more iterations — and on the smooth circuits it does the *same number
> per point* (2.00 vs 2.01), so the extra is purely the extra time points. Only `rc-pulse`
> needs more per point, 3.74 against 2.02.

## 2. Why the paper's 39% does not appear here

§4.1 reports 39% fewer time points and 17% less wall clock, and attributes it to the *lower*
bound of the two-sided band: *"the lower bound γ_min prevents step sizes from being
unnecessarily small."* That mechanism requires there to be unnecessarily small steps.

**There are not.** `benchmarks/transient_review/stage12_entry.py`: `IntegralController`'s
prediction law `h_next = h·safety·(1/err)^(1/p)` is **deadbeat** — it drives the next step's
normalised error to the fixed point `safety^p`. With `safety=0.9` and `p` recovered from the
realised step ratios as **2.00** (measured, not assumed), that fixed point is 0.81. Measured
median accepted error at `reltol=1e-6`: **0.8096 / 0.8071 / 0.8050**, IQR **0.0054 / 0.0034 /
0.0076**, with **96.4% / 94.9% / 91.2%** of accepted steps within 5% of target.

The steps already sit on the tolerance to within half a percent. **The paper is not wrong;
its baseline was different.** §4.1 records the comparison method redoing a step only above a
normalised LTE of 4.63, *without resizing toward it*. pycircuit's baseline resizes. The 39%
was available there and is not available here.

The other mechanism — rejection elimination — was measured at **0.3–3.2%** of steps, against
the stage's own 10% entry threshold, and falling as tolerance tightens. So the ceiling on
that was ~3% before any code was written.

## 3. The root cause of a month of confusion: the wrong LTE

**Equation (6), verbatim from the paper:**

    ε_m = | v_i(t_m) − v_{i,extrapolated} |

Fang's LTE is a **solution-space** quantity — computed node voltage minus a polynomial
extrapolation from previous time points, maximised over the unknowns. pycircuit computed a
divided difference of the **charge** vector, scaled by `h²` and pushed through `J⁻¹`. Same
order, different estimator. Three consequences follow, and all three mattered:

1. **Conditioning.** Charge divided differences carry repeated `1/h` factors, so as `h`
   shrinks the estimate becomes rounding noise amplified by `1/h`. An estimator that *grows*
   as the step shrinks cannot be solved for `h` — Newton runs away from the root. This is
   pinned by `test_solution_lte.py::test_it_does_not_blow_up_as_the_step_shrinks`, which
   asserts the two estimators move in opposite directions over four decades of `h`.
2. **`qᵀ` is trivial** under eq (6) — a signed unit vector on the controlling node. That is
   why §3.2 can claim *"p̄, q̄ᵀ, and d can be computed explicitly with negligible
   computational costs."* **The claim is simply false for a charge-based LTE, and that should
   have been the tell.**
3. **The controlling node moves**, which is the entire reason Figure 4 splits Newton into two
   stages. Measured here: 3 distinct nodes, switching on 0.2–13.6% of steps.

## 4. Negative results, kept

These are the ones worth having; each cost real time and each is pinned by a test so it
cannot be quietly re-adopted.

**4.1 `γ_min` is a re-parameterisation of the tolerance, and the more expensive one.**
On a deadbeat controller, raising the lower bound does not recover waste — there is none —
it moves the aim point, which is arithmetically a tolerance change. `γ_min=0.95` was paired
*in advance* with `reltol × 1.173` (the scaling that moves the aim equally):

| | steps | rejections | error |
|---|---|---|---|
| `γ_min=0.95` | 1164 (−9.6%) | **471** | 1.10× |
| `reltol × 1.173` | 1190 (−7.6%) | **20** | 1.08× |

Same trade within 2% on both axes, at 24× the rejected steps.

**4.2 The only part of the band that pays is `γ_max > 1`** — which is not the mechanism §4.1
credits. At zero accuracy cost and no step change it cuts rejections 19→1, 4→3, 48→33. That
is item 4a-bis's two-threshold controller, and it is capped by the ~3% of steps ever
rejected.

**4.3 The paper's own 0.7/3.0 band is inert here.** The aim point 0.81 lies inside it, so it
changes nothing except adding rejections (19→34, 4→25, 48→251) as steps drifting under 0.7
get redone.

**4.4 Eq (12) is catastrophically ill-conditioned for this LTE.** Eq (14) recovers `Δh` from
a denominator `qᵀdxh + d`. Those terms are the solution's sensitivity to the step size and
the extrapolation's slope — **both approximately `dv/dt`**, so their difference is the
truncation error's derivative, tiny by construction. Measured on a driven RC at `h=1.6e-7`:

    qᵀdxh = +1.818e9      d = −1.820e9      denominator = −2e6

Three digits lost, and the **sign** of the denominator decided by the cancellation. `Δh` then
saturated at the η limit with an essentially arbitrary sign, and the step drifted down four
decades while `err` sat at 0.2 — far *below* the band that should have grown it. Eq (12)
computes a small quantity as the difference of two large ones.

**This is very likely what §3.4 means by "the coupled nonlinear system sometimes is very
sensitive to the change of step size."** The shipped implementation therefore uses §3.4's
approximate Newton (eqs 17 and 18), where the step comes from the error *ratio* and no
cancellation occurs.

**4.5 The divided-difference identity is exact but buys nothing.**

    v_m − P(t_m) = v[t_m, t_{m−1}, …, t_{m−k−1}] · Π_j (t_m − t_j)

holds to machine precision, and is where the predictor and the BDF coefficients meet (both
are Newton interpolation on the same grid — for BDF, `bdf2_alphas` is the derivative at `t_m`
of the interpolant through the same nodes). It was implemented expecting better
conditioning, because the literal form subtracts two nearly-equal numbers. **It does not
help** — it relocates the cancellation into the divided-difference table, which differences
the same values. Measured on a signal shifted by a large constant, the literal form was the
*more* accurate: 9.5e-15 against 2.5e-14. Reverted; both the identity and the disproved claim
are kept as a test.

**4.6 The recorded justification for abandoning the exact update never existed.**
`_solve_coupled` carried a note that the exact (N+1) Schur update *"collapses the step
size."* Gate 12B-0 reproduced the 2026-07 code: `dh` is **identically zero in every run using
the estimator of the time**, because `E + TRTOL ≈ 0` makes `E_h ≈ 0` and the solver's
`if abs(denom) < 1e-20: dh = 0` guard fires. The code that produced that note never moved the
step at all and cannot have observed what it records. (My own first reproduction was also
invalid — it paired the stage-4 LTE formula with the pre-stage-4 charge-flavoured tolerance,
a combination that never shipped — and is withdrawn.)

## 5. Defects found, all of which produced plausible runs rather than errors

This is the recurring shape of the whole stage: **nothing raised.** Every one of these
returned a waveform.

| # | defect | how it presented |
|---|---|---|
| 1 | aim point clipped **to** `γ_min`, i.e. onto the rejection edge | 3172 rejections to accept 1187 steps, to save 7.8% |
| 2 | eq (16) damper applied to the **rejection** path | stiff ringdown crossed in 62 steps against 490, LTE reported as exactly 0 — **and the step count went DOWN**, so a step-count-only gate would have scored it as the stage's biggest win |
| 3 | step prediction used a power law where the deviation scales as `h(h+h₁)(h+h₁+h₂)` | ~`h^(k+1)` only while `h ≫ h₁`, linear when `h ≪ h₁`; mis-predicted by orders of magnitude on the shrinking side |
| 4 | predictor degree ≠ method order (Milne device) | with degree 2 against order-1 Euler, error moved **600× for a 2× step change** where the model says 4×; 2911 rejections against 3004 accepts |
| 5 | solved step discarded (`h_next` orphaned when the old loop was deleted) | **151,176 steps against 4,067**, and the tell was that the count barely moved when `reltol` changed two decades |
| 6 | no device limiting in the coupled Newton | six diode circuits returned ~0 V where the standard path gives 8.9 V — the solve "converged", to the wrong thing |
| 7 | η iterated per Newton iteration, unbounded within a time point | 0.85 per iteration over `maxiter` is seven decades; `h` reached 8.75e-15 s before giving up |
| 8 | convergence backup removed along with the LTE rejection loop | three nonlinear stress circuits failed outright |
| 9 | saturation measured on the CLAMPED step change, so a step pinned at the shrink floor reported `dh == 0.0` — indistinguishable from "stopped moving", eq (16)'s definition of converged | the first step after a pulse edge was accepted **56× too large** (2.0e-7 s where 3.55e-9 s was needed), committing a **78% single-step error**; the resulting 1.465e-2 was *identical* at reltol 1e-5 and 1e-6 |
| 10 | held steps were accepted with no error check at all — `hold_h` dropped the LTE equation *and* the test | subsumed by 9; found while chasing it |
| 11 | `dx` convergence tested against the RESIDUAL tolerance vector (`iabstol` on nodes) instead of the SOLUTION one | invisible at defaults, where `iabstol == vabstol == 1e-12`; would surface only when one is changed |
| 12 | a caller-injected step controller was silently replaced by the path's own | `tran.step_controller = X` looked honoured and did nothing — the same class as a documented feature that does not exist. Now refused with the reason |

Defect 5 is worth singling out: **it is the same defect gate 12B-0 identified in the 2026-07
code**, reintroduced by me while deleting that code. The comment in `_solve_coupled` now
records both occurrences.

Note also that fixing defect 3 alone made things *worse* (2049 → 2911 rejections). Defects 3
and 4 were independent, and the first diagnosis was incomplete.

**Defect 9 is the one worth studying**, for two reasons. First, it was invisible on every
smooth circuit: `rc-vsin` and `stiff-rlc` have no discontinuity, so nothing ever asked the
step size to fall by a factor of 56 in one time point. It took a circuit with real edges to
expose it, which is the argument for `rc-pulse` carrying a closed form.

Second, three earlier attempts at it were correct fixes to real defects and changed the
numbers *not at all* — byte-identical output across all six rows. That identity was the
signal: three different edits cannot produce the same six numbers by chance. What it meant
was that the fixed branches were not where the error came from, and only then was it worth
dumping the step sequence across the edge. The dump showed the error appearing entirely in
the *first step after* the edge and the shrink sequence pinned at exactly 0.2× — the
`MIN_SHRINK_RATIO` floor — which named the mechanism directly.

The general lesson is the one about clamps: **a clamped value must never be read as a
settled one.** Eq (16) asks whether the step size has stopped moving; a step held against a
wall has stopped moving in exactly the way that test measures, and in no other sense.
Saturation has to be measured against what the step *wants* (`h_want`, computed with the
clamps removed), not against what it got.

## 6. Design decisions that are ours, not the paper's

Recorded so they are not mistaken for transcription.

- **`τ_m`** — the paper does not define it. We reuse `LTERATIO · (reltol·ref + lte_abstol)`
  with `relref` choosing `ref`, so the coupled and standard paths are scored on one scale.
- **The extrapolation degree** follows the *active* integrator's `ORDER`, so an order drop to
  Euler drops the predictor with it. This is the Milne device; see defect 4.
- **The band target** is the geometric centre, not an edge; see defect 1.
- **`γ_min`/`γ_max` defaults are 0 and 1**, i.e. the historical one-sided test. The paper's
  0.7/3.0 are *not* adopted: they are quoted against a baseline whose redo threshold is 4.63,
  so "1.0" there is not "1.0" here. The ratio transfers; the absolute values do not.
- **A convergence backup exists.** Figure 3 forbids a backup *due to LTE*; a Newton that will
  not converge is a different failure, orthogonal to the LTE, and is retried at a smaller
  step as every production simulator does. Removing it broke three circuits.
- **`fixed_timestep` overrides the method.** The two are not really in conflict: Fang's
  method exists to *choose* the step size and `fixed_timestep` says the caller already has.
  So the grid is kept and the LTE equation is dropped on every step — the circuit is still
  solved coupled, it just has nothing left to solve for. Silently adapting anyway would
  return output points the caller did not ask for, the defect stage 4h fixed on the standard
  path.
- **`h` is floored at `minstep`** inside the coupled solve (`SchurCoupledNewton` gained an
  `hmin` argument), and the total excursion within one time point is bounded by the same
  window the standard controller allows for one step.
- **`Pulse.tr`/`tf` are floored at `MIN_EDGE = 1e-18`** so a zero-duration edge still has a
  finite derivative. Measured to cost nothing: `next_event`'s minimum-spacing guard absorbs
  the extra breakpoint, and step count, smallest step and waveform are unchanged.
- **`TLine` refuses to provide `du/dt`** rather than inheriting zero. Its source vector comes
  from history interpolated at `t − TD`, so the derivative is real and unwritten; a silent
  zero would put the coupled solve on a subtly different problem.

## 7. An invariant the method now rests on

`p` takes the right-hand limit of `du/dt` at a source's corners. That is exact rather than
approximate **only because every kink is a reported breakpoint**: the step is truncated to
land on it and `_is_first_step` re-arms, dropping the next step to backward Euler, so a step
always *starts* at a corner and moves forward into the segment ahead.

> **This was FALSE on the coupled path when first written, and is now true.** `_solve_coupled`
> had no breakpoint handling at all — no `next_event`, no truncation, no order drop — so a
> coupled step could start mid-segment and straddle a kink, which is exactly the case the
> right-hand limit assumes cannot arise. It was invisible because the only two circuits with
> closed-form references, `rc-vsin` and `stiff-rlc`, have no breakpoints.
>
> Fixing it needed two changes, not one. Truncating to the breakpoint is useless on its own,
> because `fang_timestep` *solves* for `h` and promptly replaced the truncated step: measured
> at **0 of 10 pulse edges landed on, worst miss 1.24e-7 s — the entire rise time**. A
> truncated step's size is imposed rather than free, so the LTE equation is dropped for it
> and only the circuit is solved (`hold_h`), which is the same treatment `h_clamped` gets on
> the standard path. Both paths now land on 10/10 edges to 2.7e-20 s.
> Pinned by `test_coupled_breakpoints.py`.

Because `p`'s correctness depends on that coverage being complete, it is asserted in both
directions in `test_source_derivatives.py` — every kink is reported by `next_event`, and
every reported breakpoint is a genuine jump in `du/dt`. A source that grew a new kink without
reporting it would hand `p` a derivative from the wrong segment and nothing else would
notice.

## 8. What is not done

- **Eq (12)/(13)/(14) are not used.** `SchurCoupledNewton` implements eq (12) and has the
  `hmin` floor, but the shipped path is §3.4 for the conditioning reason in 4.4. The rank-one
  LU update of §3.2 is therefore also unused.
- **`TLine` cannot be used with `coupled_lte=True`** (see §6).
- **`q̄ᵀ` and `d` are implemented, gated against finite differences, and NOT CALLED.** They
  exist for eq (12), which §3.4 replaced for the conditioning reason in 4.4; the shipped path
  needs only `p̄`, through eq (18). This is worth knowing before anyone concludes the bordered
  system is "nearly wired up" — it is, but its two cheap blocks are the only part in place and
  the expensive question (4.4) is unresolved. Found while implementing the injected-controller
  check, whose first version keyed on `lte_gradients` and would therefore have tested for a
  method nothing calls.
- **`rc-pulse` now carries a closed form** (segment-by-segment integration of the RC against
  the trapezoidal drive), so accuracy on a circuit with real edges is measurable. The first
  version of that reference was **wrong and said so loudly**: it skipped the `td` lead-in
  after the first period, where `Pulse.f` folds with `t % per` and repeats the whole shape,
  and it disagreed with a `reltol=1e-8` run by 1.0 V — full scale. Validated after the fix:
  error falls 8.1× for 9.1× more steps, i.e. `∝1/steps`, which is what a first-order method
  must give.
- **No wall-clock comparison.** Everything above counts Newton solves. §4.1's 17% is a
  runtime number and has no counterpart here.
- **The default is unchanged.** `coupled_lte` remains opt-in, and on this evidence should:
  it buys zero rejections and slightly better average accuracy for ~25% more work.
