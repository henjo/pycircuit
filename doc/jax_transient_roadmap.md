# JAX transient — limitations & refinement roadmap

**Rewritten 2026-09-01.  Every status below was checked against the code on
that date, not carried forward.**

The previous edition was a July-2026 snapshot that opened by calling this
backend a *"happy-path solver"* and listed twelve gaps, worst first.  Of
those twelve, **nine are done, two are permanently refused with cause, and
one is genuinely open** (and smaller than it was written).  It had been
stale for weeks in the expensive direction — a "not implemented" list is
what people plan work from, so a rotted one costs someone a week — and it
was made worse on 2026-09-01 by an edit that corrected item 3 and left the
other eleven, which is the sharper version of the same lesson: **correcting
one row of a status list you have just falsified is not maintenance, it is
a partial correction that certifies the rest.**

⚠ Two of the greps used to audit this document returned "present" off
**comments explaining a feature's absence**.  Grep presence is not
capability; each row below names the symbol that implements it.

## What this backend is for

`solve_batched` (P20): one compiled kernel integrating every parameter lane
of a sweep concurrently, each lane with its own adaptive — or coupled
`(x, h)` — step sequence.  The CPU deliberately has no imitation of it.
The scalar `solve()` path exists and is at parity, but it is not the reason
this backend exists, and several entries below only make sense read that
way (a `lax.cond` that is free under `jit` is not free under `vmap`).

## The 2026-07 gap list, audited

| # | gap as written | status | evidence |
|---|---|---|---|
| 1 | Device limiting (`pnjlim` inside every Newton iteration) | **DONE** | PCNR — `pcnr_inner_loop` (pair view) and `pcnr_vector_*` (device view, sec. 49), on `_pnjlim_branchless`. Also F16: the `max_dv` clamp the item complained about is now a Parameter, **default disabled** — the hardcoded 0.5 V made a 48 V rail unreachable in 100 iterations with a false "failed to converge" |
| 2 | Breakpoint order reset | **DONE** | `force_first_order` in `TransientState`, consumed by `effective_first_order` (F11/F19). Deliberately keyed on the run-global step history, not chunk-local `step_idx` — the old predicate re-dropped order at every chunk boundary, measured as chunking turning a 59-step run into 61 |
| 3 | Nonlinear continuation | **DONE** | DC: `dc_with_continuation` (junction-gmin → gshunt → Ψ-tc, P18/P25). Transient: `transient_point_rescue` (2026-09-01), gated at the dt floor, on the plain Newton **and both PCNR views**; `continuation` Parameter, default off. Source stepping stays out — see *Refused* |
| 4 | Jacobian scaling / equilibration | **REFUSED (P17)** | `nrsolver`/`scaler`/`linearsolver` are refused in `__init__`, permanently: a traced `while_loop` cannot dispatch into per-iteration Python objects. Use `Transient` |
| 5 | Order-drop on aggressive shrink | **PART DONE, part open** | see *Open* below — the zero-stability half is covered, the stalled-estimate half is not |
| 6 | Per-domain tolerances (`iabstol`/`vabstol`, `dx`+residual) | **DONE** | F6(b): `_newton_abstol`/`_newton_xtol` build the per-row vectors (iabstol on node rows, vabstol on branch rows, transposed for the update test); `conv_f` and `conv_x` are both scored, per-row, on the consistent `(F(x), dx)` pair. The old single scalar was wrong three ways at once — flavour, reference and norm |
| 7 | Integration-method selection | **DONE** | `integrator` Parameter (P6): 'gear' (default), 'euler', and 'trap' since 2026-09-01 |
| 8 | Coupled solver (Fang DAC'13) | **DONE** | `coupled_lte` Parameter, `fang_inner_loop` (P19). PCNR-inside-Fang runs too |
| 9 | Pluggable step controllers | **OPEN** | see *Open* below |
| 10 | Device bypassing (`bypass`/`bypasstol`) | **REFUSED (P13)** | A non-concept here, not a missing feature: bypassing skips evaluating quiescent elements, and a vmapped evaluation group computes all lanes of all instances in one kernel. There is nothing to skip |
| 11 | `provided_function` and `accept_step` | **DONE** | `provided_function` is supported (P11/F4 contract; must be jax-traceable, baked in at jit time). `accept_step` is not missing but *superseded*: the traced TLine ring buffer (`tline_history`/`tline_head`) replaced what the CPU's `accept_step` writes |
| 12 | DC init in `solve_batched` (was `x=0`, "for now") | **DONE** | `solve_batched` calls `dc_with_continuation` |
| — | `fixed_timestep`, `minbreak`, flexible `analysis` name | **DONE** | all three are Parameters; `fixed_timestep` gives bit-equal fixed-grid waveforms across backends |

## Genuinely open

**(A) The stalled-estimate order drop (half of old item 5).**  The item
claimed JAX "lacks Gear2's `h_curr/h_last < 0.1` protection".  Audited, that
is two different mechanisms and only one is missing:

- *Zero-stability on growth* — **covered.**  Variable-step BDF-2's parasitic
  root leaves the unit disc only on GROWTH, past
  `ZERO_STABILITY_RATIO = 1+√2 ≈ 2.414`.  The JAX controller clamps growth at
  the shared `MAX_GROWTH_RATIO = 2.0`, below the bound, and the force-accept
  path keeps the dt floor rather than growing (the CPU's own defect here was
  a force-accept path taking 10x, parasitic root 4.76).
- *The shrink branch* — **BUILT, MEASURED AND REVERTED 2026-09-01.**  Both
  CPU rules were implemented here (the ratio proxy and the `MAX_REJECT`
  cap) and both were taken out again; the numbers are in
  `benchmarks/order_drop_stalled_estimate.py` and summarised under
  *Refused* below.  In short: ported literally neither rule fires where it
  matters, and ported faithfully — to the dt floor, which is where THIS
  backend gives up — it cuts force-accepts hard and makes the answer worse
  every time.

**(B) Pluggable step controllers.**  The CPU has a `StepController` strategy
(`IntegralController`, `PIController`, `SolutionLTEController`); JAX inlines
one law — F17's `target = safety**p` band aim in `calculate_next_dt`, written
in the CPU's vocabulary so the two read as one rule.  Note this is **not**
automatically covered by P17's refusal: P17 refuses per-ITERATION Python
dispatch, whereas a step controller runs per step and its choice could be
made statically at trace time, exactly as `integrator` already is.  So this
is open, not refused — but nobody has asked for it.

**(C) Ψ-tc's continuation rung exponents — HYPOTHESIS TESTED AND DISPROVED
2026-09-01, and the real defect was elsewhere.**

The claim was that `e_start=0` / `e_max=+6` are DC-calibrated and that at the
transient dt floor the companion conductance is ~1e9, so every rung is
negligible.  Measured, that reasoning is wrong twice: the ~1e9 figure came
from a 1e-18 floor, and on the circuits where the rescue actually fires the
diagonal is ~1e-3, so the shipped exponents are too LARGE, not too small.

The measurement that matters: with gmin and gshunt disabled so Ψ-tc must
carry the rescue alone, the shipped grid rescues **1 of 4** test circuits
while a **half-decade** shift of the same exponents rescues **4 of 4**.  A
finer sweep on one circuit shows failures at offsets ≡ 0 (mod 2) — i.e.
exactly the offsets that reproduce the shipped grid of visited `g` values,
the marching step being 2 decades.  So the sensitivity is to WHICH RUNGS ARE
VISITED, not to their scale, and re-scaling would have moved the same window
somewhere else rather than fixing anything.

Chasing that turned up a genuine driver defect, now fixed — see *(F)* — but
the fix does **not** change these four results, so the grid sensitivity is
real and unexplained.  It is also **invisible in production**: gmin lands
first on every rescue demonstrated to date, so Ψ-tc never executes.  Left
alone deliberately; tuning an exponent grid to four circuits with no
mechanism would be fitting noise.

Trigger to revisit: a forced-non-converged exit that survives
`continuation=True` — i.e. a case where gmin and gshunt both fail and Ψ-tc
is actually asked to work.

**(D) Chained compact models cannot batch.**  The chained HDL path gained a
jax backend (§49.7) — 25 of 38 library classes, every real compact model,
could not be evaluated on this backend at all before.  That is a
**capability, not a speed-up**: the traced chained path is ~18x slower than
the CPU (§49.8), and §22/23's disjointness still keeps those classes out of
`solve_batched` — which is where this backend actually pays.  That
disjointness, as measured there: the models carrying real physics (PSP, the
BJTs, EKV, both MOS levels, the MESFETs) get the C backend and **cannot
batch**, while the ones that can batch (R, C, L, the ideal diode, the
passives) get no C at all — 22 / 15 / **0 both**.  (Read the two counts
carefully: 25-of-38 is today's CHAINED total from §49.7; the 22 is §22/23's
count of physics-carrying models on an earlier tree.  They answer different
questions.)  This is the largest open item on the backend by value, and it
is an `hdl.py` problem, not a transient one.

*(E) Vector PCNR inside the coupled path — CLOSED 2026-09-01, same day it
was opened.  PCNR-inside-Fang understood only the pair view; the device view,
which every real compact model produces, reached `fang_inner_loop` as the
literal string `'vector'` used as an array index.  Briefly refused at setup,
then actually wired: the Schur-reduced `(F_eff, J_eff)` is an n-sized system
whose Newton step equals predict's `dx_mna`, so fang's machinery — eq (18)'s
second solve against the same factors included — works on it unchanged, and
that argument never depended on which view produced the blocks.  Measured
against the CPU on an EKV NMOS: standard-path vector PCNR 2.000e-05, coupled
2.030e-05, i.e. as accurate as the path already trusted.*

**(F) FIXED 2026-09-01 — the ladder did not search the range it declared.**
`_adaptive_ladder_traced`'s failure path, with no rung yet landed, escalated
toward `e_max` and then refined just below `e_start`, halving until it gave
up.  The exponents it could ever try were
`{e_start, e_start+step, … e_max}` and `{e_start-1, e_start-0.5,
e_start-0.25}` — **nothing below**, however low `e_end` was set.  Measured
with a synthetic rung that converges only for `g <= 10^k`: the ladder landed
for `k >= -1` and failed for every `k <= -2`, against an `e_end` of -12
advertising a search to 1e-12.  A descent phase now walks the exponent down
to `e_end` before giving up.  It runs exactly where the old code was about to
fail, so it can only reach further; the suite was unchanged across the fix
(2,806 before and after, 2,814 once its own gates were added) and both DC
ladders use the same driver.

## Refused, with cause — do not re-propose without new measurement

- **Source stepping in the transient rescue.**  Permanent.  Scaling `u(t)`
  mid-transient would scale the integrator's companion history with it —
  ill-posed.  Ψ-tc is the mid-transient-safe substitute: its anchor is a
  state the circuit actually occupied and it scales nothing.  (DC has no
  such constraint, and the CPU DC path does source-step.)
- **P17, strategy objects** (`nrsolver`/`scaler`/`linearsolver`) — above.
- **P13, `bypass`** — above.
- **`continuation` defaulting ON.**  The CPU defaults it on because there
  the chain is Python control flow, absent from the interpreter until a
  point fails.  Here it compiles into every step of every run: measured
  ~1 s of fixed compile plus ~14% per step.  Off by default; turn it on
  when a run reports `n_forced_nonconverged > 0`.
- **The stalled-estimate order drop** (old item 5's second half).  Built
  both ways 2026-09-01, measured, reverted;
  `benchmarks/order_drop_stalled_estimate.py` reproduces it.

  Ported LITERALLY, neither CPU rule ever fires where it matters, because
  the backends give up in different places: the CPU force-accepts after
  three rejections at a point, while this loop keeps shrinking to the
  `dt_min` floor — ~130 halvings at the default 1e-18 — and force-accepts
  there.  At the floor `dt` equals the previous accepted step, so the ratio
  proxy never trips, and the point has been rejected zero times, so the
  count never reaches three.  1900 force-accepts, unchanged by either rule.

  Ported FAITHFULLY — sending a converged floor step whose LTE failed back
  for one retry at order 1 — it works as advertised and is still wrong:

  | reltol | circuit | force-accepts off → on | error off | error on | |
  |---|---|---|---|---|---|
  | 1e-6 | stiff RLC | 1900 → 516 | 5.142e-01 | 6.949e-01 | 1.35x worse |
  | 1e-6 | pulsed RC | 12 → 12 | 7.156e-04 | 1.114e-03 | 1.56x worse |
  | 1e-7 | stiff RLC | 1900 → 688 | 5.196e-01 | 6.946e-01 | 1.34x worse |
  | 1e-7 | pulsed RC | 55 → 55 | 6.935e-04 | 4.505e-03 | 6.50x worse |
  | 1e-8 | stiff RLC | 1900 → 864 | 5.208e-01 | 6.946e-01 | 1.33x worse |
  | 1e-8 | pulsed RC | 600 → 480 | 7.035e-04 | 4.498e-03 | 6.39x worse |

  Error is measured against an error-controlled reference — the same run
  with the floor lowered until nothing is force-accepted.  **Worse in 6 of
  6 configurations while cutting force-accepts by up to 3x.**

  ⚠ The transferable finding: **fewer force-accepts is not a better
  answer.**  Dropping to order 1 makes the ESTIMATE pass while making the
  actual error larger — the order-1 estimator is agreeing with a less
  accurate method about a smaller claim.  The force-accept count is a label
  on a step, not a measurement of its error, and optimising the label made
  the waveform worse every time.  The CPU benefits because its `MAX_REJECT`
  cap turns a stalled estimate into a force-accept within three retries and
  the drop prevents that; here the deep floor already prevents it, so there
  is nothing left for the drop to save and only accuracy for it to spend.

  ⚠ **A first version of this entry said the trigger was "if this backend
  ever gains a rejection cap", reasoning that a cap would produce more
  force-accepts and so give the drop something to save.  That is wrong, and
  the measurement above already refutes it**: the raised-floor runs ARE the
  force-accept-rich regime — 1900 force-accepts in 1908 accepted steps,
  more extreme than any 3-rejection cap would produce — and the drop lost
  there.  Frequency is not the variable; the drop was harmful per step it
  fired on, so firing it more often is worse, not better.

  The explanation that does fit: `h_curr/h_last < 0.1` is a DETECTOR.
  Repeated rejection is how an Integrator object, which sees only two step
  sizes, infers "we are near a discontinuity" — and that matters because a
  2nd-order LTE rests on a third divided difference, meaningless where the
  solution is not three times differentiable.  **This backend already has
  the real signal**: `force_first_order` drops the order on any step that
  lands on a declared breakpoint (F11), so on the test circuits (VPulse,
  declared edges) the corners were already integrated at order 1 before any
  of this was added.  What the rule caught instead was ordinary stiff steps
  at the dt floor — the pulsed RC has ~12 declared edges in that run and the
  rule fired 55 times — where order 2 is simply more accurate.  The proxy
  was second-guessing the estimator exactly where the estimator was right.

  So the trigger is **a discontinuity this backend has no breakpoint for** —
  a kink internal to a device model, or a piecewise I-V curve whose corner
  is in the solution but not in the time-breakpoint list.  There the proxy
  would be the only detector available, doing the job it was designed for
  rather than overriding a working estimator.

- **The coupled 'bordered' eq (12) branch.**  Measured on the CPU 2026-09-01
  before writing any JAX code (`benchmarks/coupled_bordered_gear2.py`,
  reltol 1e-5), because `bordered` was recorded as mistuned under Gear-2 and
  Gear-2 is what the JAX coupled path defaults to:

  | circuit | integrator | bordered/approx points | newton |
  |---|---|---|---|
  | pulsed | euler | 1.02x | 1.10x |
  | pulsed | **gear2** | **2.72x** | **2.24x** |
  | smooth | euler | 0.96x | 1.05x |
  | smooth | **gear2** | 0.97x | **1.18x** |

  The integrator is what mistunes it: 1.02x under Euler against 2.72x under
  Gear-2 on the same circuit, reproducing the recorded 1181-vs-350 signature.
  **And the trade is gone.**  `bordered` exists to buy fewer Newton iterations
  at the cost of a few more time points (the CPU record: 5.4% more points for
  1.6% fewer iterations).  Under Gear-2 it costs MORE on both axes, on both
  circuits.  A method that is worse in every measured dimension in the only
  configuration this backend runs is not a gap; porting it would import a
  defect and call it parity.

  What would reopen it: a re-derivation of the denominator's `w'/w` pairing
  that survives Gear-2 on the CPU.  ⚠ The obvious form of that is already a
  recorded NEGATIVE result — slicing the history to the integrator's degree
  "collapsed the bordered method to 5x the steps (3117 vs 611)" — so anyone
  picking this up should read `transient.py`'s bordered branch comments first,
  and should know that the *exact* bordered solve (the paper's `q^T dxh + d`
  denominator) is separately recorded as abandoned on evidence and must not be
  re-attempted: it is catastrophic cancellation, measured wrong-signed and 5x
  too small.

- **The rescue behind a `lax.cond`.**  The natural-looking design, and a
  trap: free under `jit`, but under `vmap` a batched predicate cannot
  branch, so `cond` becomes a `select` that evaluates both sides —
  **59x with zero lanes failing**.  Gate the ladder's loop CONDITION
  instead (0.9x, free).  Chain ladders by ANDing into the next `active`,
  never by nesting conds.

## CPU-only, with cause

**Nothing.** The list emptied on 2026-09-01. Trapezoidal and coupled +
`fixed_timestep` were ported (above); the coupled **'bordered' eq (12) branch
moved to *Refused, with cause*** on a measurement rather than a port — see
below.

*Trapezoidal integration left this list on 2026-09-01* — `integrator='trap'`.
The entry had claimed "a variable-step trap estimator has not been written";
the kernels existed and were already used by the CPU, and only the wiring was
missing. The one part that had to be written is the estimator's use of the
**charge** rather than the companion current, since differencing `g` measures
the trap recursion's undamped `(−1)^n` mode. Landing figures and the stiff-mode
ringing gate are in `doc/backend_parity_260821.md`'s P6 entry.

*coupled + `fixed_timestep` left it the same day.* `grid_locked` turned out
to be **one flag**: an over-band HELD step is normally reported unconverged so
the caller shrinks and retries, and under a caller-imposed grid shrinking is
not an option we have. Wired through `fang_inner_loop`, with the coupled
branch's `hold_h` forced true and `h_entry = grid_dt` under a fixed grid (P9's
"the grid wins" reaching the coupled path), and the accept path's `next_dt`
ordering inverted so the grid beats the "solved step carries forward" rule —
there is no solved step when the step was held. Measured: 21 steps on both
backends, waveforms agreeing to **4.4e-16** on identical points.

## Where the authoritative records live

This file is a roadmap and will rot again.  The records that are maintained
with the code are:

- **`JAXTransient`'s class docstring** — the parity ledger in force, by
  P-number.  When this file and that docstring disagree, the docstring wins.
- **`doc/backend_parity_260821.md`** — the P-entry ledger, including every
  deferral with its reopening trigger.
- **`doc/branch_review_260827.md`** — what a reviewer needs to weigh, and
  the suite totals.
- **`test_backend_conformance.py`** — the agreement, pinned rather than
  described.
