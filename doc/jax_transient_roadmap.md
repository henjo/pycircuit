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
- *The shrink branch* — **not ported.**  It has nothing to do with
  zero-stability: it is a heuristic for a STALLED ESTIMATE.  A step only
  shrinks 10x below the last accepted one after several consecutive
  rejections, and what rejects repeatedly is a 2nd-order estimate built on a
  third difference of a solution that is not three times differentiable —
  a discontinuity.  Measured value on the CPU (stiff RLC, reltol 1e-5):
  fires **3–6 times a run**, and removing it took `Gear2('ywr')` and
  `Gear2('classic')` from **0 force-accepts to 1 each**.  That is the size of
  the prize, and it is why this is worth doing rather than closing.

**(B) Pluggable step controllers.**  The CPU has a `StepController` strategy
(`IntegralController`, `PIController`, `SolutionLTEController`); JAX inlines
one law — F17's `target = safety**p` band aim in `calculate_next_dt`, written
in the CPU's vocabulary so the two read as one rule.  Note this is **not**
automatically covered by P17's refusal: P17 refuses per-ITERATION Python
dispatch, whereas a step controller runs per step and its choice could be
made statically at trace time, exactly as `integrator` already is.  So this
is open, not refused — but nobody has asked for it.

**(C) Ψ-tc's continuation rung exponents.**  `e_start=0`, `e_max=+6` are
calibrated for DC diagonals.  At the transient dt floor the companion
conductance is ~1e9, so every rung of that ladder is negligible against it.
gmin and gshunt carry **every** rescue demonstrated so far, so this is a
third resort that has not been reached.  Trigger: a forced-non-converged
exit that survives `continuation=True` on the PCNR path.

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
