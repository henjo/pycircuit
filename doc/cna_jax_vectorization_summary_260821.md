# Branch summary: `cna-jax-vectorization` — 2026-08-21

A merge-readiness overview for the repo owner.  The branch is complete by
its own roadmap: every scheduled item is executed, gated, and recorded;
the suite stands at **1005 passed, 3 skipped, 3 xfailed**; and the
feature the branch exists for now carries a measured payoff number.
This document is the map — each section points at the detailed record.

## What the branch delivers

**The JAX backend with `solve_batched`** — vmapped, jit-compiled
transient analysis running every lane of a parameter sweep concurrently,
each lane with its own swept parameters and its own DC bias.  Measured
payoff (`doc/batched_sweep_260821.md`): on a 512-lane rectifier corner
sweep, **22.5×** over the honest CPU loop (54.7 s vs 2.4 s warm), with
the crossover at ~16 lanes and the JAX wall nearly flat in N — the
speedup keeps growing with lane count at this circuit size.  Below ~16
lanes the CPU loop wins and remains the honest recommendation.

**CPU/JAX feature parity, tested rather than claimed**
(`doc/backend_parity_260821.md`, P1–P25, all executed).  Both backends
now share: parameter vocabulary and defaults, Gear-2 as the shipped
default (P22's state-row mask made the coupled path safe at order 2),
integrator selection, all three `relref` modes with the unit-group
split, the LTE acceptance band, `fixed_timestep`, `provided_function`,
real initial conditions (the dead `node.ic` path deleted), TLine with
quadratic monotone-limited interpolation, the Fang coupled estimator and
PCNR (both traced), and the voltage/current excursion checks
(`max_dv_step`/`max_di_step` with `'auto'` from sampling theory).
Deliberate divergences are documented where users look, each with its
reason (P13 bypass is a non-concept under vmap; P17 strategy objects
cannot survive a traced loop — permanent, stated).

**A complete DC/transient continuation family** (P18 + P25, both
backends): junction-gmin (the physical homotopy) → gshunt (the
singular-G rescue) → source stepping (CPU DC only) → pseudo-transient
Ψtc (the basin-respecting last resort), all driven by one shared
adaptive exponent-space ladder (halving-on-failure, escalation, and an
exactly-zero final rung — only a pure converged solution is ever
accepted).  The CPU transient's failed-point rescue uses the same chain
minus source stepping (which would scale companion history).  SER
δ-adaptation was prototyped, measured, and rejected with the numbers
recorded; the JAX transient-point rescue is deferred with its full
architectural reasoning and a concrete reopening trigger (both in the
P25 entry).

**The repair campaign that preceded parity**
(`doc/transient_review_260820.md`, F1–F19 + hygiene, all executed;
baselines and exit gate in `doc/phase0_baseline_260821.md`): the
`solve_batched` crash and silent-override defects, the JAX Newton
criterion (wrong flavour, wrong reference, wrong norm), the missing
breakpoint order drop, the t=0 result point, the safety factor, the
0.5 V clamp, and the estimator/order coupling — each with the test that
failed before the fix.  Two standing instruments guard the defect
classes the review identified: the dead-knob scan (dead arguments) and
the cross-backend conformance harness (behavioral divergence), both in
the suite.

## Numbers a reviewer might want first

- rc conformance bound at matched order: **4.3e-6** (Phase-3 exit gate).
- Coupled+Gear-2 rectifier: 259 points vs Euler's 769 at equal accuracy.
- JAX PCNR cost: **+29 % wall** vs the CPU's +60–80 %/iteration (the
  figure did not transfer to vmapped execution, as suspected).
- TLine cross-backend: matched line **bit-close (5.6e-16)**; the JAX
  DC-stamp defect it exposed (24.5 V on a matched 1 V line, untested
  since landing) is fixed and gated.
- Batched sweep: crossover ~16 lanes, 22.5× at 512, ≤1 mV cross-backend
  final-sample agreement at matched order.

## What is deliberately NOT here (each with its trigger, recorded)

- JAX transient-point continuation rescue (P25 entry: the ladder would
  have to be traced into the step body; trigger = a forced-non-converged
  chunk exit on a circuit the CPU rescue completes).
- Batch compaction for finished lanes (invisible in the flat benchmark
  wall; trigger = many-lane profiles showing dispatch cost).
- `Sin.next_event` quarter-period breakpoint coalescing (trigger = a
  profile showing breakpoint cost on sinusoidal drives).
- `bordered` Fang and curvature equidistribution beyond `max_dv_step`
  (each "until separately justified").

## Review order suggestion

1. `doc/backend_parity_260821.md` — the parity contract and its
   execution record (the closest thing to a change log by intent).
2. `doc/transient_review_260820.md` — why the repairs were made, with
   measurements and the meta-review record.
3. `doc/phase0_baseline_260821.md` + `doc/batched_sweep_260821.md` —
   the before/after yardsticks.
4. The suite: `pytest pycircuit/circuit/tests/` (1005 green as of this
   summary; the conformance harness and dead-knob scan are the standing
   gates).

Merging is the owner's call; nothing in this branch assumes it.
