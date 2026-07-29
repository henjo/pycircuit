# Symbolic distortion on decision diagrams — implementation plan

**Status: planned, not started.**

**Goal, stated by the maintainer: run distortion analysis on bigger circuits.**
That is the objective the stages are ordered against; reducing `sympy.cancel`
cost (§8.3) is the presenting symptom, not the target. The two are not the
same problem — see `distortion_ddd_conclusions.md` §8, which measures that the
dominant cost *changes regime* with circuit size.

Reasoning and feasibility measurements behind this plan:
`doc/distortion_ddd_conclusions.md`. Read that first — it contains the
measurement that decides whether the plan is worth running at all, and the one
trap that would make a prototype look hopeless for the wrong reason.

Answers `doc/distortion_mimo_plan.md` §8.3, which established that symbolic
expression *size* is polynomial in truncation order and established **nothing
about cost** — while `sympy.cancel` failed to finish in 900 s at `U^7`.

## 1. Scope, and what stays untouched

**Nothing in this plan changes the existing implementation.** New module,
new tests, new documents, everything under a `distortion_ddd` label:

```
pycircuit/circuit/distortion_ddd.py          new
pycircuit/circuit/tests/test_distortion_ddd.py   new
doc/distortion_ddd_conclusions.md            written
doc/distortion_ddd_plan.md                   this file
doc/src/circuit/distortion_ddd.rst           new, if stage D is reached
```

`distortion.py` is **read-only for this work**. The reason is not tidiness:
the existing path is the oracle. Refactoring it to share code with its own
checker would destroy the independence that makes the cross-check worth
anything.

**Out of scope, with reasons:**

- **Speeding up the numeric path.** Even the fast DDD route costs ~3 us per
  solve against numpy's sub-microsecond dense solve. The numeric path stays
  numpy. **Reconsider if** a circuit appears where the matrix is large and
  sparse enough that a shared-cofactor solve wins — not before measuring it.
- **Replacing `graded_response_mimo`.** The graded ring, the consistency
  filter and the composition bookkeeping are correct and gated. This work
  changes *what a coefficient is made of*, not how the recurrence runs.
- **Two-tone, initially.** Single tone first; the index machinery is shared
  and two-tone adds nothing to the representation question. **Reconsider if**
  stage C lands cheaply.

## 2. The gates, and the unusual luxury behind them

Every stage below is checked against the **existing implementation**, which is
gated on four published closed forms (eqs. 41, 42, 48, 52), an independent ODE
integration (magnitude, phase to 0.0000 deg, waveform to 2e-10) and a 2-D
numerical Fourier extraction.

This is a better oracle than the DDD work itself had. Its calibration boundary
is explicit: for named op-amps "a published `|DDD|` pins an order of magnitude,
not a number". Here a symbolic result can be required to agree to floating
point with something already known to be right.

**So no gate in this plan is a self-comparison, and none is a plausibility
argument.**

## 3. Stages

### Stage A — the cheap alternative, measured first

**Not a DDD stage.** Represent each graded coefficient as a sum of products of
*unevaluated* factors — no `cancel`, no `expand`, no premultiplication — and
evaluate through a single memo, exactly the discipline that gave 1000x in the
feasibility spike.

**Gate.** On the two-node symbolic system that defeated `cancel`: produce the
`U^7` third-harmonic coefficient and evaluate it, in **under 60 s**, agreeing
with the current numeric path to 1e-12.

**This stage can end the plan.** If a plain factored ring makes `U^7`
tractable, DDD is an optimisation of a solved problem. **Record that as the
answer to §8.3 and stop** — it is a result, not a failure, and it is the
outcome that costs least to discover.

### Stage B — a DDD-backed solve

Replace the `solve` callable with one built on `DDDFamily` plus `eval_roots`:
build the cofactor table once, memoise per distinct harmonic frequency,
recombine per right-hand side.

**Gate.** `graded_response_mimo` driven through it must reproduce the existing
biquad and amplifier results **to floating point** — the same eq. (48) and
eqs. (41)/(42) comparisons, at the same tolerances. Not "close": the
arithmetic is the same arithmetic, so anything above ~1e-12 means the
representation is wrong.

**Guard against the trap that would fake a failure here:** `eval` in a loop
costs 1000x `eval_roots` with a memo. If stage B looks hopeless, check that
first.

### Stage C — sharing across harmonics via `s_expand`

Use `s_expand` so one construction serves every harmonic, rather than one
family per frequency.

**Gate.** Same numerical agreement as stage B, plus a **measured** reduction in
distinct vertices against stage B's one-family-per-frequency arrangement. If
sharing does not reduce the count, say so: it would mean the harmonics have
less structure in common than the shared pencil suggests, which is worth
knowing and contradicts the motivating argument.

### Stage C2 — the lever for bigger circuits, before any new diagram

**Ahead of the exotic options deliberately.** `HierarchicalDDD` and
`suppression_order` already exist and are already validated on a circuit far
larger than anything the distortion work has run: 1040 → 156 vertices on the
µA741, and 11 088 → 26 vertices with 86 s → 0.48 s on noise. Suppressing
internal nodes is a property of **the matrix**, and the distortion path solves
the same matrix at every harmonic — so this should transfer with no new
representation at all.

**Gate.** Distortion on a circuit substantially larger than the biquad — take
one from `benchmark_circuits.py` — with the harmonics reproduced against the
existing numeric path to 1e-10, and a **measured** vertex reduction from
suppression. If suppression does not reduce the count here although it does
for the plain determinant, that is a finding about the graded solve worth
recording, not a failure to hide.

**This stage is the one that most directly serves the stated goal.** If it
works, bigger circuits become reachable without inventing anything.

### Stage D — is it worth it?

**Gate, stated in advance so it is not rationalised afterwards.** The DDD path
is worth keeping if it either

1. makes a symbolic truncation order tractable that the stage A
   representation cannot reach at all, or
2. reduces representation size by more than 3x on a circuit larger than the
   biquad, at equal accuracy.

If it does neither, **record it as a negative result and do not ship it as a
default.** The DDD work kept five such results deliberately; this plan should
be equally willing to produce a sixth.

Scale matters for this stage: the 2x2 biquad has five shared vertices and
proves nothing. Use the existing `benchmark_circuits.py`.

### Stage E — a *graded* multiroot diagram, if C2 leaves a gap

The natural extension of what exists rather than a new family:
`SExpandedDDD` keys its memo on `(rows, cols, power-of-s)`; the distortion
work grades by `(harmonic, power-of-drive)`. One root per graded key, memo key
`(rows, cols, harmonic, power)`, is the same multiroot construction applied to
the expansion variable this analysis already carries.

**Gate.** Fewer distinct vertices than stage C at equal accuracy, on the
stage C2 circuit. **Do not start this before C2 reports** — if hierarchy
already reaches the target circuit sizes, this is effort spent on a solved
problem.

Polynomial-shaped diagrams (TED, `*BMD`) are deliberately *not* a stage. They
are the right shape for the numerator polynomial and are discussed in
`distortion_ddd_conclusions.md` §9.2, but they would be new machinery, their
compactness under a variable set that grows with truncation order is
unmeasured, and **no number should be taken from that literature without the
paper in hand.**

## 4. What would make this fail, stated in advance

- **The products, not the determinants, are the bottleneck.** After `p`
  passes a coefficient is a sum of products of up to `p` cofactor/determinant
  factors at different frequencies. DDD compacts and shares each factor; it
  does not stop that tree growing. If the growth is dominated by the product
  structure rather than by the factors, stages B and C will show correct
  results and no win.
- **Too few distinct frequencies to share across.** Sharing pays when many
  harmonics draw on one construction. A `U^3` truncation touches three.
- **Float terminals.** The DDD work's exact-rational concerns do not apply
  (negative result P3: a diagram never multiplies entries symbolically), but
  conditioning at high harmonic index is unmeasured here.

## 5. Outcomes

_To be filled in as stages complete. Negative results recorded here, not
edited out — the DDD work kept five and they are among its more useful
output._

| Stage | Outcome |
|---|---|
| A — factored ring, no DDD | not started |
| B — DDD-backed solve | not started |
| C — `s_expand` across harmonics | not started |
| C2 — hierarchy, for bigger circuits | not started |
| D — worth it? | not started |
| E — a graded multiroot diagram | not started, and optional |

## 6. Carried-over facts that bear on this

- **`eval_roots`, never `eval`, in a loop** — 3.07 us against 3.05 ms,
  measured. The single most important implementation note here.
- **`DDDFamily` uses `'min-degree'`, not `'auto'`**, because it addresses
  cofactors by *original* row and column and `'auto'` permutes the matrix.
  Do not change it.
- **Complex `s` is fine** — 2.5e-15 against numpy at three harmonics.
- **Keep factors unmultiplied as tuples**, as `DDDCombination` does: a
  premultiplied product is a sympy expression and costs a full substitution.
- **`sympy.cancel` does not finish in 900 s at `U^7`** on a two-node system.
  That is the number stage A has to beat, and the reason this plan exists.
